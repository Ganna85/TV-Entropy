function nbins = get_N_bins(Gamma)
            nbins = 1;
            for t = 1:size(Gamma,2)-1
                if ~isequal(Gamma(:,t),Gamma(:,t+1))
                    nbins = nbins + 1;
                end
            end
        end