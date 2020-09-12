% BPSK UW Acoustic Communication Equalization script
% Author: Pedro Córdoba González
%

close all; clear all;
addpath('Simulated Channel Response'); % We add to the path the folder with the CRs
addpath('Functions'); % We add to the path the folder with the CRs
%% Parameters
%Communication Parameters
M = 2;          % Order of the modulation
L_data = 100;  % Number of transmitted symbols
L_lea = 20;    % Number of known symbolsh
Fs_sym = 500;   % Symbol Frequency
SNR = 0:2:20;   % Signal to Noise Ratio
Packets = 100; % Number of transmitted symbols

%Channel Parameters
Fs_h=1e4;       % Sample frequency of Channel Impulse Response
Fs_c=3e4;       % Sample Frequency of Chirp
CRfile='Frequency_Response_sim_dir_45-55kHz_25Hz_60s_0.05s_395_5_25_OK.mat';
Channel_data=load(CRfile); % Data simulated with Stojanovic script
Lf=401; Lt_tot=3603; T_SS=60; T_tot=3*T_SS;
fmin=5e3; % minimum frequency [Hz]
B=10e3; % bandwidth [Hz]
df=25; % frequency resolution [Hz], f_vec=fmin:df:fmax;
dt=50e-3; % time resolution [seconds]
T_SS=60; % coherence time of the small-scale variations [seconds]
shift=10; skip=10;

% Equalizer Parameters
nTaps = 20; %Length of the equalizer
delay = 0;
%DFE
dfe = comm.DecisionFeedbackEqualizer('Algorithm','LMS', ...
    'NumForwardTaps',20,'NumFeedbackTaps',10,'StepSize',0.03,'ReferenceTap', 1);

%We create the progress bar and initiate the variables
f = waitbar(0,'Calculating BERs...');
ber_nf = zeros(Packets, length(SNR));
ber_zf = zeros(Packets, length(SNR));
ber_mmse = zeros(Packets, length(SNR));
ber_dfe = zeros(Packets, length(SNR));
k_ind = zeros(Packets, length(SNR));

%we create the chirp to estimate the channel
t = 0:1/Fs_c:1-1/Fs_c;
swept = chirp(t,0,t(end),B)';

%We start the simulation
for i=1:Packets
    %% Channel adjustment
    k = round(rand()*(length(Channel_data.hmat)-1))+1; %Number of CR selected
    k_ind(i) = k; %We save the index of the channel used for tests
    hmat = Channel_data.hmat;
    H = Channel_data.H_LS;
    h_raw = circshift(hmat(:, k), shift); % From all the CR we select a random one
    h_raw = h_raw/norm(h_raw); %Normalization of the CR
    
    
    %Resample of the Channel's Response from the first arrival
    [p,q] = rat(Fs_sym / Fs_h);
    [m,ind] = max(abs(h_raw(1:50)));        %calculation of the first arrival
    % h_sym = resample(h_raw(ind:end),p,q);
    h_sym = h_raw(ind:q:end);
    Lsym = length(h_sym);
    
    %Adjustment for Channel Estimation - We resample the Channel's Response to fit the chirp
    %We resample the Channel's Response to fit the chirp
    [p,q] = rat(Fs_c / Fs_h);
    H_raw = H(:, k);
    H_raw_res = resample(H_raw,p,q);
    h_raw_res = circshift(ifft(H_raw_res), shift);
    % h_raw_res = ifft(H_raw_res);
    [m,ind] = max(abs(h_raw_res(1:200)));        %calculation of the first arrival
    h_raw_res = h_raw_res(ind:end);
    h_raw_res = h_raw_res/norm(h_raw_res); %Normalization of the CR
    
    
    %% Communication
    % Generation of the random training data
    data = randi([0 1],L_data,1);
    
    % bpsk mapper
    data_mod = pskmod(data, M);
    
    % Calculation of the symbols received
    data_r_nonoise = conv(h_sym,data_mod);
    
    % We conv the two signals
    swept_r = cconv(h_raw_res,swept);
    
    for j=1:length(SNR)
        %% Channel Estimation
        % We are going to use a chirp to estimate the channel response as we would
        % do in the real communication
        
        %We add the different noises to the signal received
        swept_r = awgn(swept_r,SNR(j));
        
        %We calculate the estimated Channel's Response
        Y = swept_r(1:Fs_c);
        Sxy = conj(swept).*Y;
        Sxx = conj(swept).*swept;
        Syy = conj(Y).*Y;
        H_r = Sxy./Sxx;
        h_r = ifft(H_r);
        
        %We resampled to the symbol Frequency
        [p,q] = rat(Fs_sym / Fs_c);
        h_sym_r = h_r(1:q:end);
        h_sym_r = h_sym_r(1:Lsym);
        
        
        %% Data reception
        % Adding noise to the signal affected by the channel and
        % calculating variance to get the MMSE filter
        data_r = awgn(data_r_nonoise,SNR(j));
        var_s = var(data_r_nonoise);
        var_n = 10^((10*log10(var_s)-SNR(j))/10);
        
        % Demodulation of the data without filtering
        data_demod = pskdemod(data_r, M);
        
        [~, ber] = biterr(data_demod(1:L_data), data);
        ber_nf(i,j)=ber;
        
        % Calculating ZF filter
        H = fft(h_sym_r);
        h_eq = ifft(1./H_r);
        h_eq = h_eq(1:q:end);
        h_eq =  h_eq(1:Lsym);
        data_eq = filter(h_eq,1,data_r);
        %Demodulation of ZF filtered signal
        data_demod_eq = pskdemod(data_eq, M);
        [~, ber] = biterr(data_demod_eq(1:L_data), data);
        ber_zf(i,j)=ber;
        
        % Calculating MMSE filter
        h_eq = ifft(conj(H_r)./((abs(H_r).^2)+(var_n/var_s)));
        %         h_eq = ifft(conj(H)./((conj(H).*H)));
        h_eq = h_eq(1:q:end);
        h_eq =  h_eq(1:Lsym);
        data_eq = filter(h_eq,1,data_r);
%         data_eq = conv(h_eq,data_r);
        % Demodulation of MMSE filtered signal
        data_demod_eq = pskdemod(data_eq, M);
        [~, ber] = biterr(data_demod_eq(1:L_data), data);
        ber_mmse(i,j)=ber;
        
        %Filtro DFE
        data_eq = dfe(data_r, data_mod);
        % Demodulation of DFE filtered signal
        data_demod_eq = pskdemod(data_eq, M);
        [~, ber] = biterr(data_demod_eq(1:L_data), data);
        ber_dfe(i,j) = ber;
    end
    waitbar(i/Packets,f,'Calculating BERs...');
end

close(f)

%We average all the BERs calculated
BER = mean(ber_nf);
BER_zf = mean(ber_zf);
BER_mmse = mean(ber_mmse);
BER_dfe = mean(ber_dfe);

mBER = median(ber_nf);
mBER_zf = median(ber_zf);
mBER_mmse = median(ber_mmse);
mBER_dfe = median(ber_dfe);

figure
plot(SNR, BER,'-*')
hold on
plot(SNR,BER_zf,'-*')
plot(SNR,BER_mmse,'-*')
plot(SNR, BER_dfe,'-*')
legend('No filter','ZF','MMSE','DFE')
title('Averaged Bit Error Rate (BER)','fontsize', 16);
xlabel('SNR (dB)', 'fontsize', 12), ylabel('Rate', 'fontsize', 12)

figure
plot(SNR, mBER,'-*')
hold on
plot(SNR,mBER_zf,'-*')
plot(SNR,mBER_mmse,'-*')
plot(SNR, mBER_dfe,'-*')
legend('No filter','ZF','MMSE','DFE')
title('Median Bit Error Rate (BER)','fontsize', 16);
xlabel('SNR (dB)', 'fontsize', 12), ylabel('Rate', 'fontsize', 12)

figure
histogram(ber_mmse(:,7),10)
title(['Histogram MMSE BER for SNR= ', num2str(SNR(7))],'fontsize', 14);
xlabel('BER', 'fontsize', 12), ylabel('Number of cases', 'fontsize', 12)