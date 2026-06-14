function [r1, noise] = delay_time_transmission(N, M, sigma_2, s, gs, L_set)

        r1=zeros(N*M,1);
        noise= sqrt(sigma_2/2)*(randn(size(s)) + 1i*randn(size(s)));
        for q=1:N*M
            for l=(L_set+1)
                if(q>=l)
                    r1(q)=r1(q)+gs(l,q)*s(q-l+1);
                end
            end
        end
        r1=r1+noise;
end