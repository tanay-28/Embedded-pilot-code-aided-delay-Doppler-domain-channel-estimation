function [X_tilda, s] = otfs_modulation(X, Fn, N, M)
        X_tilda=X*Fn';% delay-time transmit symbols
        s = reshape(X_tilda,N*M,1);
end
