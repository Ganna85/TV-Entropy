function [hidden] = get_Hidden(K, T, LP, x, val);
               
    Gamma_hat=zeros(K, LP.tBins);
    for k=1:K
        Gamma_hat(k,:)=x((k-1)*(2*LP.tBins-1)+1:k*(2*LP.tBins-1)-LP.tBins+1);
    end
    
    %Gamma out of "Gamma_hat"
    gamma=zeros(K,T);
    for t=1:LP.tBins
        gamma(:, LP.pos_switch(t)+1:LP.pos_switch(t+1)) = Gamma_hat(:,t)*ones(1,LP.pos_switch(t+1) - LP.pos_switch(t));
    end
    
    hidden.gamma = gamma;
    hidden.nbins = get_N_bins(gamma);
    hidden.x = x;
    hidden.val = val;
    hidden.LP = LP;
               
    