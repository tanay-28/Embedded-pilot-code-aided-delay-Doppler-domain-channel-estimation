clearvars;
%% OTFS parameters%%%%%%%%%%
% N: number of symbols in time
N = 16;
% M: number of subcarriers in frequency
M = 64;
% M_mod: size of QAM constellation
M_mod = 4;
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));

%% delay-Doppler grid symbol placement
% ZP length  should be set to  greater than or equal to maximum value of delay_taps
length_ZP = M/16;
M_data = M-length_ZP;
data_grid=zeros(M,N);
data_grid(1:M_data,1:N)=1;
Nsym = sum(sum(data_grid));

%% SNR 
SNR_dB = 0:2:30;
SNR = 10.^(SNR_dB/10);
sigma_2 = (abs(eng_sqrt)^2)./SNR;
N_fram = 10;
%% Normalized DFT matrix
Fn=dftmtx(N);  % Generate the DFT matrix
Fn=Fn./norm(Fn);  % normalize the DFT matrix
Fm=dftmtx(M);  % Generate the DFT matrix
Fm=Fm./norm(Fm);  % normalize the DFT matrix
%% generating bits for data frame
data_sample=qammod(reshape(randi([0,1],Nsym*M_bits,1),M_bits,Nsym), M_mod,'gray','InputType','bit');
X_sample = Generate_2D_data_grid(N,M,data_sample,data_grid);
s_sample = reshape(X_sample*Fn',N*M,1);
sample_avg = mean(abs(s_sample).^2);
%%
l_p = M/2;
l_tau = length_ZP-1;
k_p = N/2;
k_v = 2;

code_length = 17;
c1 = exp(1j*p3(code_length,1));

half_code_length = (code_length-1)/2;
l_gaurd = half_code_length;

% Pilot symbol constant values
Xp1 = sqrt(code_length);
Xp2 = 1;

% initializing error count variables
errors_LMMSE = zeros(N_fram,length(SNR_dB));
errors_LMMSE_b = zeros(N_fram,length(SNR_dB));
errors_LMMSE_i = zeros(N_fram,length(SNR_dB));

% 3GPP synthetic channel configuration
taps = 4;% no. of channel taps

% range of absolute value of channel coefficients
lower_threshold = 1;
higher_threshold = 1.4;
%% channel gen
% uncomment for same channel for each frame
% [G, gs, H_tf, H_dd, H_dt, H, H_tilda, P, delay_taps, Doppler_taps, L_set] = Synthetic_channel_gen(M,N,Fn,Fm,taps, l_tau, k_v, lower_threshold, higher_threshold);
%%
for iesn0 = 1:length(SNR_dB)
    for ifram = 1:N_fram
%% channel gen
% different channel for each frame
[G, gs, L_set, H_tf, H_dd, H_dt, H, H_tilda(:,:,ifram,iesn0), P, delay_taps(:,ifram,iesn0), Doppler_taps(:,ifram,iesn0), chan_coef_2(:,ifram,iesn0)] = Synthetic_channel_gen(M,N,Fn,Fm,taps, l_tau, k_v, lower_threshold, higher_threshold,"simulate");

%% setting the data grid
l_tau = max(max(delay_taps(:,ifram,iesn0)));
k_v = round(max(max(abs(Doppler_taps(:,ifram,iesn0)))));
data_grid(l_p-l_tau-half_code_length:l_p+l_tau+half_code_length,k_p-2*k_v:k_p+2*k_v) = zeros(2*l_tau+code_length,4*k_v +1);
N_syms_perfram = sum(sum(data_grid));
N_bits_perfram1 = N_syms_perfram*M_bits;

data_grid2=data_grid;

N_syms_perfram2 = sum(sum(data_grid2));
N_bits_perfram2 = N_syms_perfram2*M_bits;

%% uncoded transmit data gen

[X, data, trans_info_bit] = transmit_data_gen(N, M, N_bits_perfram1, N_syms_perfram, M_bits, M_mod, data_grid);
 X = X/sqrt(2);% setting the abs value of the symbols to 1

%% coded transmit data gen

[X2, data2, trans_info_bit2] = transmit_data_gen(N, M, N_bits_perfram2, N_syms_perfram2, M_bits, M_mod, data_grid2);
 X2 = X2/sqrt(2); 

%% adding pilot symbols

 X(l_p,k_p)=Xp1;
 X2(l_p-half_code_length:l_p+half_code_length,k_p) = c1*(Xp2);

%% OTFS modulation

[X_tilda, s] = otfs_modulation(X, Fn, N, M);
[X_tilda_2, s2] = otfs_modulation(X2, Fn, N, M);

 %% delay-time domain transmission

[r1, noise] = delay_time_transmission(N, M, sigma_2(iesn0), s, gs, L_set);
[r2, noise2] = delay_time_transmission(N, M, sigma_2(iesn0), s2, gs, L_set);

%% OTFS demodulation

[Y_tilda, Y, Y_tf] = otfs_demodulation(r1, N, M, Fn, Fm);
[Y_tilda2, Y21, Y_tf2] = otfs_demodulation(r2, N, M, Fn, Fm);

%% coded channel estimation

threshold_1 = code_length-2;
threshold_2 = 1.4;
[Y1(:,:,ifram,iesn0),Y_corr, est_delay_tap, est_Doppler_tap, tap_value, tap_value_compensated, est_L_set_2, est_taps, taps_in_row(iesn0,:)] = barker_coded_pilot_estimation(Y21,l_p,k_p,l_tau,l_gaurd, k_v, Xp2, c1, code_length, half_code_length, threshold_1, threshold_2, phase);


    end
end

%% plotting PAPR curves
papr2 = mean(papr2,1);
papr1 = mean(papr1,1);
% 
[N1,X3] = histcounts(papr1, 100,'Normalization','count');
[N2,X4] = histcounts(papr2, 100,'Normalization','count');
figure
semilogy(X3(1:end-1),1-cumsum(N1)/max(cumsum(N1)),'-s','LineWidth',2)
hold on
semilogy(X4(1:end-1),1-cumsum(N2)/max(cumsum(N2)),'-x','LineWidth',2)
grid on
% xlim([3 5])
xlabel('PAPR')
ylabel('CCDF')
lgd2 = legend('uncoded pilot symbol','barker coded pilot symbol');
lgd2.Location = 'southwest';
lgd2.FontSize = 14;

%% plotting the recovered channel
figure;
% subplot 211
% image(abs(est_H_dd),'CDataMapping','scaled');
% subplot 212
% image(abs(H_dd),'CDataMapping','scaled');
bar3(abs(H_dd))
xlabel('Doppler')
ylabel('delay')
zlabel('Amplitude')
set(gca,'XTickLabel',-k_v:k_v)

%% plotting the cross-correlation of the output
figure;
bar3(abs(Y_corr))
xlabel('Doppler')
ylabel('ACR')
zlabel('Amplitude')
set(gca,'XTickLabel', -k_v:k_v)
set(gca,'YTickLabel')

%% plotting the BER curves
figure;
semilogy(SNR_dB,mean(errors_LMMSE),'-s','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,mean(errors_LMMSE_b),'-x','LineWidth',2,'MarkerSize',8)
hold on
semilogy(SNR_dB,mean(errors_LMMSE_i),'-x','LineWidth',2,'MarkerSize',8)
grid on
xlabel('SNR (dB)')
ylabel('BER')
lgd1 = legend('uncoded pilot symbol',' barker coded pilot symbol','ideal channel estimation');
lgd1.Location = 'southwest';
lgd1.FontSize = 14;
%%
save('otfs_channel_data.mat', 'delay_taps', 'Doppler_taps', 'chan_coef_2','Y1');