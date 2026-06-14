function [G, gs, L_set, H_tf, H_dd, H_dt, H, H_tilda, P, delay_taps, Doppler_taps, chan_coef_2] = Synthetic_channel_gen(M,N,Fn,Fm,taps, l_tau, k_v, lower_threshold, higher_threshold, mode)

        if mode == "test1"
               
        chan_coef_2 = ones(1,taps);
        delay_taps = [0 1 2 3];
        Doppler_taps = [-1 0 1 2];

        elseif mode == "test2"

        chan_coef_2 = ones(1,taps);

        delay_taps = randperm(l_tau+1, taps)-1;
        Doppler_taps = randperm(2*k_v,taps)-k_v;

        elseif mode == "test3"

        chan_coef_2 = ones(1,taps);

        delay_taps = randperm(l_tau+1, taps)-1;
        Doppler_taps = [-2 0 1 2];
            
        elseif mode == "test4"

        chan_coef_2 = ones(1,taps);

        delay_taps = [0 2 1 3];
        Doppler_taps = randperm(2*k_v,taps)-k_v;

        elseif mode == "simulate"

        chan_coef=1/sqrt(2)*(randn(1,20*taps)+1i.*randn(1,20*taps));
        chan_coef_2 = chan_coef(find((abs(chan_coef)>lower_threshold)&(higher_threshold>abs(chan_coef)),taps));

        delay_taps = randperm(l_tau+1, taps)-1;
        Doppler_taps = randperm(2*k_v,taps)-k_v;

        end

        L_set=unique(delay_taps);
        [G,gs]=Gen_time_domain_channel(N,M,taps,delay_taps,Doppler_taps,chan_coef_2);
        [H_tf]=Generate_time_frequency_channel_ZP(N,M,gs,L_set);
        [H,H_tilda,P]= Gen_DD_and_DT_channel_matrices(N,M,G,Fn);
        H_tf = (H_tf).';
        H_dd = fftshift(Fm'*H_tf*Fn);
        H_dt = (Fm'*H_tf).';
        H_tf = (H_tf).';
        end