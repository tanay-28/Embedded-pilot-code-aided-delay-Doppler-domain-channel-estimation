function [hcap,hcap_compensated,h_est_dd,channel_coef_compensated,est_delay_taps,est_Doppler_taps,est_L_set,est_taps] = embedded_pilot_channel_estimation(M,N,Y,l_p,k_p,k_v,l_tau,sigma_2,Xp, phase)
%% extraction of the channel image
        Y1 = Y(l_p:l_p+l_tau,k_p-k_v:k_p+k_v);
        Tau = 3*sqrt(sigma_2);
        b = zeros(l_tau+1,2*k_v+1);
        scaling = 2^(0.5*log2(M)+2);  
%% Thresholding        
        for k=1:1 + 2*k_v
            for l=1:l_tau+1
                      if abs(Y1(l,k))>Tau
                         hcap(l,k) = Y1(l,k)./(Xp);% delay-Doppler channel estimate
                         b(l,k) = 1;
                      else
                         hcap(l,k) = 0;
                      end
            end
        end
 %% Doppler compensation
   d = exp(1i*phase);
   for it1=1:2*k_v+1
       hcap_compensated(:,it1) = hcap(:,it1)./(d^(it1-1-k_v));
   end
%% delay-Doppler reconstructed channel
   h_est_dd = zeros(M,N);
   h_est_dd(l_p+1:l_p+l_tau+1,k_p-k_v+1:k_p+k_v+1) = hcap*scaling;
%% 
   channel_coef_compensated=[];
   for it21=1:l_tau+1
   channel_coef_compensated= cat(2,channel_coef_compensated,nonzeros(fliplr(hcap_compensated(it21,:))).');
   end

est_Doppler_taps = [];
est_delay_taps = [];
for u=1:l_tau+1
    elements = sort(find(abs(hcap(u,:))>0),'descend');
    est_Doppler_taps = cat(2,est_Doppler_taps,elements-k_v-1);
    for f = 1:length(elements)
    est_delay_taps = cat(2,est_delay_taps,u-1);
    end
end

est_L_set = unique(est_delay_taps);
est_taps = length(est_delay_taps);
end