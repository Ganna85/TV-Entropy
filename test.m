%% run TV-Entropy model on S&P500 sample used in the paper

%% load SPX data
SPX = load('/data/SPX/x.mat');

%% get scaled [-1 1] sample data
xt = SPX.xt';

%% Estimate 2-regime model
%% set the numebr of regimes
K = 2;
%% set the number of switches
C = [2];
%% set the L1-bounds fore each regime
L1_bounds = [Inf Inf];
%% set the num of moments for the ME problem
k = [6 6];
%% set the number of the anneiling steps
anneil = 2;

[result] = main(xt, K, C, L1_bounds, k, anneil, 40, 1e-06, 1e-06, 'SPX_test');

gamma = zeros(1,length(xt));
%% plot densities of the regimes against Gaussian
for i=1:result.BestModel.K
    regime_data = xt(result.BestModel.hidden.gamma(i,:) == 1);
    vars(i) = var(regime_data);
    figure;
    histogram(regime_data, 'Normalization', 'pdf', 'DisplayName', ['regime ', num2str(i)]);
    hold on;
    scatter(xt, exp(-result.BestModel.resid(i,:)), 'r', 'DisplayName', 'MaxEnt');
    norm_dens = normpdf(xt,mean(regime_data),std(regime_data));
    hold on;
    scatter(xt, norm_dens, 'k', 'DisplayName', 'Gaussian');
    title(['Regime ',num2str(i)]);
    
    ind = find(result.BestModel.hidden.gamma(i,:)==1);
    gamma(ind) = i;
end

%% plot regime-switching path
figure;
plot(xt, 'DisplayName', 'SPX');
hold on;
plot(gamma, 'r', 'Linewidth', 3, 'DisplayName', 'switching path');

