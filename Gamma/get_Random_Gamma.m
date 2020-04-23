function [hidden] = get_Random_Gamma(K, T, C, LP);
    %generate random gamma respecting constraints in LP
    %MIX-integer liner programming problem with gurobi solver    
    if (K == 1)
        hidden.gamma = ones(1,T);
    else
        max_trys = 500;
        for i = 1:max_trys
            temp = -rand(K, T);
            c_Theta = get_C_Theta(LP, K, temp);
            
            [x, val, exitflag, result] = linprog_gurobi(c_Theta, LP.Aneq, ...
                LP.bneq, LP.Aeq,LP.beq,zeros(K*(2*LP.tBins-1),1),...
                ones(K*(2*LP.tBins-1),1), [], 0);

            %sometimes gurobi returns continuous solution, not OK!
            if length(x(x(x > 0) < 1)) > 0
                display('gamma-step: get_Random_Gamma gamma is not binary');
            end
                [hidden] = get_Hidden(K, T, LP, x, val);

                if ~all(sum(hidden.gamma,2))
                    display('gamma-step: get_Random_Gamma - empty cluster');
                else
                    break;
                end
            end
    end