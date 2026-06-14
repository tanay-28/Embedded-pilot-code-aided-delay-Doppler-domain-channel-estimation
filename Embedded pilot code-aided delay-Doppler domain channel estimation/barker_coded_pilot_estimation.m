function [Y1,Y_corr, est_delay_tap, est_Doppler_tap, tap_value, tap_value_compensated, est_L_set, est_taps, taps_in_row] = barker_coded_pilot_estimation(Y,l_p,k_p,l_tau,l_gaurd, k_v,Xp, barker_seq, barker_length, half_barker_length, threshold_1, threshold_2, phase)

 Y1 = Y(l_p-l_gaurd:l_p+l_gaurd+l_tau,k_p-k_v:k_p+k_v).';
 b1 = zeros(2*l_gaurd+l_tau+1,2*k_v+1);
 taps_in_row = zeros(2*barker_length+l_tau-1,1);
 j=1;
 k=1;
 f=1;
%% cross correlation
        Y_corr = xcorr2(Y1,barker_seq).';

        for itc=1:size(Y_corr,1)
            Y2 = Y_corr(itc,:)/Xp;
            [max_values_Y2, max_index_Y2] = maxk(Y2,2);% extracting the two largest values
            max_value_Y2 = max(max_values_Y2);
            
         if(abs(max_value_Y2)>threshold_1)% checking if the maximum value is above a threshold (delay contains a tap)
         
             
             if abs(max_values_Y2(1))/abs(max_values_Y2(2))>threshold_2 % checking whether there is one tap
                taps_in_row(itc) = 1;
                tap_value(k) = max_value_Y2/barker_length;
                ind1 = find(Y2==max_value_Y2);
                est_Doppler_tap(f) = ind1-k_v-1;
                idx(itc) = ind1;
                b1(itc-l_gaurd,ind1)=1; 
                k=k+1;
                f=f+1;
             else                       % checking whether there are two taps
                taps_in_row(itc) = 2;
                [temp,ind] = maxk(Y2,2);

                est_Doppler_tap(f)=ind(1)-k_v-1;
                est_Doppler_tap(f+1)=ind(2)-k_v-1;

                tap_larger = Y2(ind(1))/barker_length;%larger tap
                tap_smaller = Y2(ind(2))/barker_length;%smaller tap

                if(abs(tap_larger-Y1(ind(1),itc-half_barker_length))>abs(tap_smaller-Y1(ind(2),itc-half_barker_length))) %
                    tap_value(k) = tap_larger;
                    tap_value(k+1) = tap_smaller;
                else
                    tap_value(k) = tap_smaller;
                    tap_value(k+1) = tap_larger;
                end

                b1(itc-l_gaurd,ind(1))=1;
                b1(itc-l_gaurd,ind(2))=1;
                k=k+2;
                f=f+2;
            end
            
         end
    end
   j=1;
%% Calculating delay taps
est_delay_tap =[];
   for itc1=1:size(b1,1)
       for itd1 = 1:size(b1,2)
           if b1(itc1,itd1)==1
               est_delay_tap(j) = itc1-l_gaurd-1;
               j=j+1;
           end
       end
   end

   est_taps = sum(taps_in_row);
   %% Doppler compensation
   d = exp(1i*phase);
   for it2 = 1:est_taps
       tap_value_compensated(it2)=tap_value(it2)./(d^est_Doppler_tap(it2));
   end
% est_L_set_2 = [];
est_L_set = unique(est_delay_tap);
end