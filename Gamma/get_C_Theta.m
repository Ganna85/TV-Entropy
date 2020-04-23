function C_Theta = get_C_Theta(LP, K, Resid);
    %Residuals as matrix "G_hat"
    G_hat = zeros(K,LP.tBins);
    for j=1:LP.tBins
        G_hat(:,j)=sum(Resid(:,LP.pos_switch(j)+1:LP.pos_switch(j+1)),2);
    end
               
    %Residuals as vector "C_hat"
    C_Theta=zeros((2*LP.tBins-1)*K,1);
    for k=1:K
        C_Theta((k-1)*(2*LP.tBins-1)+1:k*(2*LP.tBins-1)-LP.tBins+1,1) = G_hat(k,:);
    end
    