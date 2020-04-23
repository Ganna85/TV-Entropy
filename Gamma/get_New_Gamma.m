function [hidden, result] = get_New_Gamma(K, T, C, resid, hidden);
    %generate random gamma respecting constraints in LP
    %MIX-integer liner programming problem with gurobi solver    
    if (K == 1)
        hidden.gamma = ones(1,T);
        result = [];
        result.runtime = 0;
        result.status = '';
        hidden.val = 0;
    else
        c_Theta = get_C_Theta(hidden.LP, K, resid);
        
        [x, val, exitflag, result] = linprog_gurobi(c_Theta, hidden.LP.Aneq, hidden.LP.bneq,...
            hidden.LP.Aeq, hidden.LP.beq, zeros(K*(2*hidden.LP.tBins-1),1),ones(K*(2*hidden.LP.tBins-1),1),...
            hidden.x, hidden.val);
                          
                      
        if length(x(x(x > 0) < 1)) > 0
            display('gamma-step: gamma is not binary');
            x = round(x);
        else
           %x = round(x);
           [hidden] = get_Hidden(K, T, hidden.LP, x, val);
           if ~all(sum(hidden.gamma,2))
                display('gamma-step: empty cluster');
           end
        end
    end