function [x0_norm, Mom, LB, x0, X_, Aineq, bineq] = get_inputs(eps, gamma, X, k, bGARCH);
    %get indexes of points active in this regime
    points_i = find(gamma);
    
    %calcualte k moments
    Mom = sum(X(points_i,1:k),1)/sum(gamma);
    X_ = X(points_i,1:k);
    data = X(points_i,1);
    
    %Calculate lagrange multipers for normal distribution,
    %to use as starting points for estimation
    if (bGARCH)
        m = 0;
        v = 1;
    else
        m = mean(data);
        v = var(data);
    end
    
    L0 = -log(1/sqrt(2*pi*v)) + m^2/(2*v);
    L1 = -m/v;
    L2 = 1/(2*v);
    
    if (~isinf(eps))
        %start values for LM - close to normal
        x0_norm = zeros(2*k,1);
        x0_norm(1:2) = [L1 L2]';

        %Boundary positive slack vars
        %LB = [-Inf*ones(k-1,1); 0; zeros(k,1)];

        %last LM and all slack vars >=0
        LB = [-Inf*ones(k-1,1); zeros(k+1,1)];
        
        %x start
        x0 = zeros(2*k,1);
        
        %Ax<=b
        pI = eye(k);
        nI = -1*eye(k);
        Aineq = [pI nI; ...
        nI(1:end-1,:) nI(1:end-1,:); ...
        zeros(1,k) ones(1,k)];
        bineq = [zeros(2*k-1, 1); eps];
                
    else
        %start values for LM - close to normal
        x0_norm = zeros(k,1);
        x0_norm(1:2) = [L1 L2]';

        %Boundary positive slack vars
        %LB = [-Inf*ones(k-1,1); 0; zeros(k,1)];

        %last LM >=0
        %LB = [-Inf*ones(k-1,1); 0];
        LB = [-Inf*ones(k,1)];
        
        %x start
        x0 = zeros(k,1);
        
        Aineq = [];
        bineq = [];
                
    end
    
end