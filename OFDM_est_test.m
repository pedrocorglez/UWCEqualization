% BPSK UW Acoustic Communication Equalization script
% Author: Pedro Córdoba González
%

close all; clear all;
addpath('Simulated Channel Response'); % We add to the path the folder with the CRs
addpath('Functions'); % We add to the path the folder with the CRs
%% Parameters
%Communication Parameters
L_data = 100; % Number of transmitted symbols
L_sym= 128;
L_lea = 1000;
Fs_sym = 250;   % Symbol Frequency
SNR = 0:2:20;   % Signal to Noise Ratio

%Channel Parameters
Fs_h=1e4;       % Sample frequency of Channel Impulse Response
Fs_c=3e4;       % Sample Frequency of Chirp
CRfile='Frequency_Response_sim_seq_45-55kHz_25Hz_60s_0.05s_395_5_25_OK.mat';
Channel_data=load(CRfile); % Data simulated with Stojanovic script
Lf=401; Lt_tot=3603; T_SS=60; T_tot=3*T_SS;
fmin=45e3; % minimum frequency [Hz]
B=10e3; % bandwidth [Hz]
df=25; % frequency resolution [Hz], f_vec=fmin:df:fmax;
dt=50e-3; % time resolution [seconds]
T_SS=60; % coherence time of the small-scale variations [seconds]
shift=10; skip=10;

%Modulation Parameters
M = 2; % Modulation order for QPSK
phase = 0;
nfft  = 128;
cplen = 16;
B_mod = 6e3;

% OFMD Parameters
K = 128; %number of OFDM subcarriers

CP = K/4; %length of the cyclic prefix: 25% of the block

P = 17; %number of pilot carriers per OFDM block


allCarriers = 1:K; % indices of all subcarriers ([1, 1, ... K])

pilotCarriers = 1:(K/(P-1)):K; %Pilots is every (K/P)th carrier.
%For convenience of channel estimation, let's make the last carriers also be a pilot
pilotCarriers = [pilotCarriers, K];

% data carriers are all remaining carriers
dataCarriers = allCarriers;
dataCarriers(pilotCarriers)=[];
dataL_sym = length(dataCarriers); % number of payload bits per OFDM symbol

%We create the progress bar and initiate the variables
f = waitbar(0,'Calculating BERs...');
ber_nf = zeros(L_data, length(SNR));
ber_zf = zeros(L_data, length(SNR));
ber_mmse = zeros(L_data, length(SNR));
% ber_dfe = zeros(L_data, length(SNR));
OFDM_symbol = zeros(Lf,1);

for i = 1:L_data
    % Channel acquisition
    k = round(rand()*(length(Channel_data.hmat)-1))+1; %Number of CR selected
    H_LS = Channel_data.H_LS;
    hmat = Channel_data.hmat;
    H_raw = H_LS(:,k);

    % Channel adjustment
    H_carriers = ((length(H_raw)-1)/2)+1-K/2:((length(H_raw)-1)/2)+K/2;
    % H = H_raw(H_carriers);
    % h_raw = circshift(ifft(H), shift);   % From all the CR we select a random one
    % h = ifft(H);
    % % h_raw = h_raw/norm(h_raw);            % Normalization of the CR
    h_raw = circshift(hmat(:, k), shift); % From all the CR we select a random one
    [m,ind] = max(abs(h_raw(1:25)));        % Calculation of the first arrival
    h = [h_raw(ind:end); zeros(ind-1,1)];
    H=fft(h);
    H=H(H_carriers);
    h=ifft(H);
    
    %     % Channel adjustment
%     H_carriers = ((length(H_raw)-1)/2)+1-K/2:((length(H_raw)-1)/2)+K/2;
%     % H = H_raw(H_carriers);
%     % h_raw = circshift(ifft(H), shift);   % From all the CR we select a random one
%     % h = ifft(H);
%     % % h_raw = h_raw/norm(h_raw);            % Normalization of the CR
%     h_raw = circshift(hmat(:, k), shift); % From all the CR we select a random one
%     [m,ind] = max(abs(h_raw(1:25)));        % Calculation of the first arrival
%     h = [h_raw(ind:end); zeros(ind-1,1)];
%     H=fft(h);
% %     H=H_raw;
% %     H = zeros(length(H_raw),1)
% %     H(H_carriers)=H_raw2(H_carriers);
% %     h=ifft(H);
    
    %% Communication
    % We create the bits for each symbol
    data = randi([0 M-1],dataL_sym,1); %K random data
    
    % We modulate the data in PSK symbols
    data_psk = pskmod(data, M, phase);
    pilot_psk = pskmod(zeros(P,1),M, phase); % The known value each pilot transmits
    
    symbol(dataCarriers) = data_psk;
    symbol(pilotCarriers) = pilot_psk;
%     OFDM_symbol(H_carriers) = symbol;
    
    % We calculate the symbol in time
%     symbol_t = ifft(OFDM_symbol);
    symbol_t = ifft(symbol);
    
    % Channel convolution
    symbol_t_r = conv(symbol_t, h);
%     symbol_t_r = symbol_t_r(1:Lf);
    symbol_t_r = symbol_t_r(1:K);
    
    for j=1:length(SNR)
        
        % Adding white noise
        symbol_t_r_noise = awgn(symbol_t_r, SNR(j), 'measured');
        
        % Back to Frequency domain
%         OFDM_symbol_r = fft(symbol_t_r_noise);
%         symbol_r = OFDM_symbol_r(H_carriers);
        symbol_r = fft(symbol_t_r_noise);
        
        % Channel estimation
        pilots_r = symbol_r(pilotCarriers); %Extraction of the pilots from the received symbol
        H_est_pilots = pilots_r ./ pilot_psk; % divide by the transmitted pilot values
        
        x = 1:P;
        xq = 1:(P-1)/K:P-(P-1)/K;
        % %interpolación lineal
        % H_est_a = interp1(x, abs(H_est_pilots), xq, 'linear');
        % H_est_p = interp1(x, angle(H_est_pilots), xq, 'linear');
        % H_est = H_est_a .* exp(1i*H_est_p);
        %interpolación cuadrática
%         H_est_a = interp1(x, abs(H_est_pilots), xq, 'spline');
%         H_est_p = interp1(x, angle(H_est_pilots), xq, 'spline');
%         H_est = H_est_a .* exp(1i*H_est_p);
        %interpolación cuadrática
        H_est_a = interp1(x, abs(H_est_pilots), xq, 'cubic');
        H_est_p = interp1(x, angle(H_est_pilots), xq, 'cubic');
        H_est = H_est_a .* exp(1i*H_est_p);
        
        % Equalization zFE
        symbol_zfe = symbol_r ./ H;
        
        %Equalization MMSE
        symbol_mmse = (conj(H).*symbol_r) ./ abs(H).^2;
        
        % Demodulation without equalization
        data_r = pskdemod(symbol_r(dataCarriers), M, phase);
        [~, ber_nf(i,j)] = biterr(data_r, data);
        
        % Demodulation with equalization
        data_r_zfe = pskdemod(symbol_zfe(dataCarriers), M, phase);
        [~, ber_zf(i,j)] = biterr(data_r_zfe, data);
        
        % Demodulation with equalization
        data_r_mmse = pskdemod(symbol_mmse(dataCarriers), M, phase);
        [~, ber_mmse(i,j)] = biterr(data_r_mmse, data);
        
    end
    waitbar(i/L_data,f,'Calculating BERs...'); 
end

close(f)

%We average all the BERs calculated
BER = mean(ber_nf);
BER_zf = mean(ber_zf);
BER_mmse = mean(ber_mmse);
% BER_dfe = mean(ber_dfe);

mBER = median(ber_nf);
mBER_zf = median(ber_zf);
mBER_mmse = median(ber_mmse);
% mBER_dfe = median(ber_dfe);

figure
plot(SNR, BER,'-*')
hold on
plot(SNR,BER_zf,'-*')
plot(SNR,BER_mmse,'-*')
% plot(SNR, BER_dfe,'-*')
legend('No filter','ZF','MMSE','DFE')
title('Averaged Bit Error Rate (BER)','fontsize', 16);
xlabel('SNR (dB)', 'fontsize', 12), ylabel('Rate', 'fontsize', 12) 

figure
plot(SNR, mBER,'-*')
hold on
plot(SNR,mBER_zf,'-*')
plot(SNR,mBER_mmse,'-*')
% plot(SNR, mBER_dfe,'-*')
legend('No filter','ZF','MMSE','DFE')
title('Median Bit Error Rate (BER)','fontsize', 16);
xlabel('SNR (dB)', 'fontsize', 12), ylabel('Rate', 'fontsize', 12) 

figure
histogram(ber_mmse(:,7))
title(['Histogram MMSE BER for SNR= ', num2str(SNR(7))],'fontsize', 14);
xlabel('BER', 'fontsize', 12), ylabel('Number of cases', 'fontsize', 12) 

figure(4)
plot(allCarriers, real(H(H_carriers)),'g')
hold on
grid on
stem(pilotCarriers, real(H_est_pilots))
plot(allCarriers, real(H_est),'b')
legend('Correct','Pilots','Linear int.')
xlim([0 K+1])