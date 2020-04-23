function [LP] = getLPConstraints(K, T, Ct, num_fe ); 
    %builds matrices representing constraints in gamma-optimization
    LP = [];
    LP.tBins = num_fe - 1;
    LP.pos_switch = round(0:(T)/(LP.tBins):T);
    %==============================================================
            
    % Matrix for the inequalities
    B = zeros(LP.tBins-1, LP.tBins);
    B(:,1:LP.tBins-1) = -eye(LP.tBins-1);
    B(:,2:LP.tBins) = B(:,2:LP.tBins) + eye(LP.tBins-1);
    C = -eye(LP.tBins-1);
    D = -B;
            
    E = zeros(1,LP.tBins);
    F = ones(1,LP.tBins-1);
            
    Aneqi = [B C ; D C ;E F]; % block i-th of Aneq
            
    helpMeCell = cell(1,K);
    [helpMeCell{:}] = deal(Aneqi);
    LP.Aneq =  blkdiag(helpMeCell{:});
    LP.bneq = repmat([zeros(2*LP.tBins-2,1); Ct], K,1);
            
    % Matrix for the equalities
    LP.Aeq = repmat([eye(LP.tBins) zeros(LP.tBins,LP.tBins-1)],1,K);
    LP.beq = ones(LP.tBins,1);