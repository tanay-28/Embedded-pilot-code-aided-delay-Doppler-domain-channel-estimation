function [sn_block_est, Gn, Rn, X_est_dd] = generate_Gn(N,M,noise_var,r,gs,L_set,Fn)

for n=1:N    
    for m=1:M
        for l=L_set+1
            if(m>=l)
                Gn(m,m-l+1)=gs(l,m+(n-1)*M);
            end
        end
    end
    rn=r((n-1)*M+1:n*M);    
    Rn=Gn'*Gn;
    sn_block_est(:,n)=(Rn+noise_var.*eye(M))^(-1)*(Gn'*rn);
end

X_est_dd = sn_block_est*Fn;
end