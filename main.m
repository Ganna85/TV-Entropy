function [result] = main(xt, K_range, C_range, eps_range, k, anneil, subsp_iter, eps_subsp, tresh_signific, name);
%% xt         - input data [1,T];
%% K_range    - # of clusters range [1,K], e.g., [1 2 3];
%% C_range    - # of switches range, e.g., [10:10:T]
%% eps_range  - L1 regularization parameters range, e.g., [1e-06, 1e-05, 1e-04, 1e-03]
%% k          - # of moments for Max Ent estimation, e.g., 6
%% anneil     - # of anneiling steps, e.g., 10
%% subsp_iter - # of subspace iterations e.g., 40
%% eps_subsp  - treshold of subspace, e.g., 1e-06
%% tresh_signific - treshold for choosing significant parameters in Lasso, e.g., 1e-06
%% name - name of the output file
%% ------------------------------------
%% result.Input       - Input parameters
%% result.Candidates  - Best results per combination accross anneilings [1 x #combinations]
%% result.BestModel   - best BIC model
%%
%% where  #combinations = #K_range*#C_range*#eps_range;


%% log the outputs into a text file
%diary(['Results/', name]);

result = [];
inp.xt = xt;
inp.K_range = K_range;
inp.C_range = C_range;
inp.eps_range = eps_range;
inp.k = k;
inp.anneil = anneil;
inp.subsp_iter = subsp_iter;
inp.eps_subsp = eps_subsp;
inp.tresh_signific = tresh_signific;

result.Inputs = inp;
result.Candidates = cell(0);
result.BestModel = cell(0);

%% save seed
reset(RandStream.getGlobalStream);
rng('shuffle');
result.Seed = rng;

%% check input data
checkinputs(inp);

%% Create the power matrix
X = [];
X = inp.xt';
for i=2:max(inp.k)
    X = [X inp.xt'.^i]; 
end

ind_candidates = 1;
T = length(xt);

%% set options for the solver
options = optimset('Display','off', 'MaxFunEvals', 5000, 'MaxIter', 5000 , 'TolFun', 1e-12,...
                                'TolX', 1e-10, 'TolCon',1e-12, 'GradObj','on', 'GradConstr','on', 'Algorithm','interior-point', 'DerivativeCheck','off');
                            

%% estimate models for all paramer combinations
for K = K_range(1:end)
        for C = C_range(1:end)
            display('-----------------------------------------------------------------------------------------------------------------');
            display(['####### Parameters combination: K = ', num2str(K), ', k = [', num2str(k),'], C = ', num2str(C), ', eps = [', num2str(eps_range),']']);
            
            %% Best result (w.r.t. BIC) for this parameter combination will be stored here
            Best = cell(0);
            
            %% Build constraints matrices for LP problem
            LP = get_LP_Constraints(K, T, C, T);
            ind_anneil = 1;
            
            %% idea of anneiling that we try to estimate model with same
            %% parameters "anneil" number of times from different initial
            %% values of gamma, and then to choose the best one
            for an = [1:anneil]
                %% flag to check if subspace converged
                bConverged = false;
                
                %% STEP0: generate initial random gamma (using gurobi solver for Mix-Integer problem)
                hidden = [];
                hidden = get_Random_Gamma(K, T, C, LP);
                
                try
                %% run subspace (STEP1 and STEP2) until the convergence OR until the number of
                %% iterations reached subsp_iter
                    for s = 1:subsp_iter 
                        optRes = [];
                        optRes.acf = Inf;
                        optRes.hidden = hidden;
                        flag = [];
                        
                        %% STEP1: (lambda-step)
                        tstart = tic;
                        temp = 0;
                        for i=1:K
                            param = [];
                            
                            %% get input values for optimization
                            [x0_norm, param.Mom, LB, x0, ~, Aineq, bineq] = get_inputs(eps_range(i), optRes.hidden.gamma(i,:), X, k(i), false);
                            something_wrong = false;
                            
                            %% run lambda-optimization
                            try
                                [y, ~, param.fval, param.flag, param.output, f, param.g, param.h, param.time] = ...
                                    MLE_test(k(i), param.Mom, options, x0_norm, x0, LB, Aineq, bineq, eps_range(i), -1, 1);
                            catch
                                display('Error: lambda-optimization error');
                            end

                            param.Lagr = y(1:k(i));
                            param.Lagr_long = y;
                            
                            %% find normalization constant
                            param.Z = integral(@(x) exp(-x.^(1:k(i))*param.Lagr), -1, 1, 'ArrayValued', true, 'AbsTol', 1e-12, 'RelTol', 1e-9);

                            %% find L0
                            param.L0 = log(param.Z);

                            %% find the norm (pdf should integrate to 1)
                            param.norm = integral(@(x) exp(-x.^(0:k(i))*[param.L0; param.Lagr]), -1, 1, 'ArrayValued', true, 'AbsTol', 1e-12, 'RelTol', 1e-9);

                            %% check the moments error, should be 0 for L1 = INF.
                            param.I = [];
                            for ind=1:k(i)
                                param.I(ind) = integral(@(y) y.^ind*exp(-y.^(1:k(i))*param.Lagr), -1, 1, 'ArrayValued', true, 'AbsTol', 1e-12, 'RelTol', 1e-9);
                                param.I(ind) = param.I(ind)/param.Z;
                            end
                            param.err = param.I - param.Mom;

                            %% find value of L1-norm
                            param.RP = sum(abs(param.Lagr));
                            param.eps = eps_range(i);
                            param.k = k(i);
                            
                            flag = [flag param.flag];
                            %% store Max Ent density parameters for this regime
                            optRes.params{i} = param;
                            %% compute the residium
                            optRes.resid(i,:) = (X(:,1:k(i))*param.Lagr + ones(T,1)*param.L0)';
                            temp = temp + param.fval;
                        end
                        %% this measures the time
                        time = toc(tstart);
                        
                        %% compute value of the loss-function L
                        optRes.acf = sum(sum(optRes.resid .* optRes.hidden.gamma,1) ,2);
                        display(['param optim: acf ', num2str(optRes.acf), ', time = ', num2str(time), ' exitflag ', num2str(flag), ', objval = ', num2str(temp,10)]);

                        if (s == 1)
                           prev = Inf;
                        end

                        %% see if acf converged or if acf grows
                        if (prev < optRes.acf)
                            display(['Error: Growing ACF after parameter estimation detected, returning previous best result ', num2str(prev)]);
                        else
                            if (abs(prev - optRes.acf) < eps_subsp)
                                %% converged! - > break
                                bConverged = true;
                                break;
                            end
                       end

                       %% STEP2: (gamma-step)
                       %% using gurobi
                       [hidden, res] = get_New_Gamma(K, T, C, optRes.resid, optRes.hidden);
                       
                       %% Calculate new value of L
                       acf = sum(sum(optRes.resid .* hidden.gamma,1) ,2);  
                       display(['gamma optim: acf ', num2str(acf), ', time = ', num2str(res.runtime), ' exitflag ', res.status, ', objval = ', num2str(hidden.val,10)]);

                       %% check again if acf grows
                       if optRes.acf < acf
                            display(['Error: Growing ACF after gamma estimation detected, returning previous best result ', num2str(optRes.acf)]);
                            break;
                       end

                       optRes.acf = acf;
                       optRes.hidden = hidden;
                       prev = optRes.acf;
                    end
                catch ex
                % Get last segment of the error message identifier.
                    fprintf(['>> annealing step ', num2str(an),...
                        ', subspace step ',int2str(s),...
                        ', invalid solver state detected (', ...
                          ex.message,'). Iteration deleted.\n']);
                end
                            
                if (bConverged)
                    %% calculate Information Criterion for this model
                    optRes.K = K;
                    optRes.C = C;
                    optRes.T = T;
                    optRes.k = k;
                    optRes.eps_range = eps_range;
                    optRes = compute_IC(optRes, tresh_signific);
                    %% save result in a list
                    Best{ind_anneil} = optRes;
                    display(['Anneiling step #', num2str(an), ' converged after ', num2str(s),' iterations with acf = ', num2str(optRes.acf)]);
                    %% and start new anneiling: doing the same thing for the same parameters
                    %% but with different initial gamma
                    ind_anneil = ind_anneil + 1;
                else
                    %% if did not converge after subsp_iter iterations
                    display(['Anneiling step #', num2str(an), ' didnt converge after ', num2str(s),' iterations with acf = ', num2str(optRes.acf), '. Iteration deleted']);
                end
            end
            
            %% in the list of BEST, that contains results for all anneiling steps, we look for model with best(min) BIC,
            min_BIC = Best{1}.ic.BIC;
            for j = 1:length(Best)
                if (Best{j}.ic.BIC <= min_BIC)
                    min_BIC = Best{j}.ic.BIC;
                    min_BIC_j = j;
                end
            end
            
            result.Candidates{1,ind_candidates} = Best{min_BIC_j};
            display(['Best BIC result is acf:', num2str(Best{min_BIC_j}.acf), ' and BIC:', num2str(Best{min_BIC_j}.ic.BIC)]);
            ind_candidates = ind_candidates + 1;    
        end
end

%% find best(min) BIC among all combinations or parameters
display('--------------------------------------------------------------------------------------------------');
minBIC = result.Candidates{1}.ic.BIC;
BIC_j = 1;
for j = 2:length(result.Candidates)
    if (result.Candidates{j}.ic.BIC < minBIC)
        minBIC = result.Candidates{j}.ic.BIC;
        BIC_j = j;
    end
end

result.BestModel = result.Candidates{BIC_j};

display(['Best BIC model: K = ', num2str(result.BestModel.K), ', C = ', num2str(result.BestModel.C),...
    ', k = [', num2str(result.BestModel.k), '], eps = [', num2str(result.BestModel.eps_range), '], BIC = ',...
    num2str(result.BestModel.ic.BIC)]);

%% save the output file
% save(['Results/', name], 'result');
% diary OFF;


function [optRes] =  compute_IC(optRes, tresh_signific);
    T = optRes.T;
    K = optRes.K;
            
    logLike = -optRes.acf;
            
    %% compute level of persitence
    hidParNum = 0;
    if (K > 1)
        hidParNum = optRes.hidden.nbins - 1;
    end
    
    %% compute the numner of sign parameters in ME densitites
    modParNum = sum(optRes.k);
    for i=1:K
         for (j = optRes.k(i):-1:1)
            if (abs(optRes.params{i}.Lagr(j)) <= tresh_signific)
                modParNum = modParNum - 1;
            else
                break;
            end
         end
    end
    
    optRes.ic.modParNum = modParNum;
    optRes.ic.hidParNum = hidParNum;

    %% estimate all number of params
    number_param = hidParNum + modParNum;
            
    %% before computing IC we check if overfitting occures in one of
    %% the clusters
    time_in_cluster = sum(optRes.hidden.gamma,2);
    if all(time_in_cluster > modParNum/K)
        %% compute BIC
        optRes.ic.BIC = -2*logLike + log(T) * number_param;
        optRes.ic.true_BIC = -2*logLike + log(T) * (modParNum + 1);
    else
        disp('Overfitting, force BIC = inf')
        optRes.ic.BIC = Inf;
        optRes.ic.true_BIC = Inf;
    end      