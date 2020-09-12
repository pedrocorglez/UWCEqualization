% QPSK UW Acoustic Communication Equalization script
% Author: Pedro Córdoba González
%

close all; clear all;
addpath('Simulated Channel Response'); % We add to the path the folder with the CRs
addpath('Functions'); % We add to the path the folder with the CRs
%% Parameters
%Communication Parameters
M = 4;          % Order of the modulation
L_data = 100;  % Number of transmitted symbols
L_lea = 20;    % Number of known symbolsh
Fs_sym = 500;   % Symbol Frequency
SNR = 0:2:20;   % Signal to Noise Ratio
Packets = 1000; % Number of transmitted symbols

%Channel Parameters
Fs_h=1e4;       % Sample frequency of Channel Impulse Response
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
EQ = 2;
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

%We start the simulation
for i=1:Packets
    %% Channel adjustment
    k = round(rand()*(length(Channel_data.hmat)-1))+1; %Number of CR selected
    hmat = Channel_data.hmat;
    h_raw = circshift(hmat(:, k), shift); % From all the CR we select a random one
    h_raw = h_raw/norm(h_raw); %Normalization of the CR
    
    
    %Resample of the Channel's Response from the first arrival
    [p,q] = rat(Fs_sym / Fs_h);
    [m,ind] = max(abs(h_raw(1:50)));        %calculation of the first arrival
    % h_sym = resample(h_raw(ind:end),p,q);
    h_sym = h_raw(ind:q:end);
    Lsym = length(h_sym);
    
    %% Communication
    % Generation of the random training data
    data = randi([0 3],L_data,1);
    
    % bpsk mapper
    data_mod = pskmod(data, M, pi/4);
    
    % Calculation of the symbols received
    data_r_nonoise = conv(h_sym,data_mod);
    
    for j=1:length(SNR)
        % Adding noise to the signal affected by the channel and
        % calculating variance to get the MMSE filter
        data_r = awgn(data_r_nonoise,SNR(j));
%         var_n = 10^(-SNR(j)/10);
%         var_s = var(data_r_nonoise);
        var_s = var(data_r_nonoise);
        var_n = 10^((10*log10(var_s)-SNR(j))/10);
        
        % Demodulation of the data without filtering
        data_demod = pskdemod(data_r, M, pi/4);
        [~, ber] = biterr(data_demod(1:L_data), data);
        ber_nf(i,j)=ber;
        
        % Calculating ZF filter
        H = fft(h_sym);
        h_eq = ifft(1./H);
        h_raw = h_eq/norm(h_eq); %Normalization of the CR
        data_eq = conv(h_eq,data_r);  
        %Demodulation of ZF filtered signal
        data_demod_eq = pskdemod(data_eq, M, pi/4);
        [~, ber] = biterr(data_demod_eq(1:L_data), data);
        ber_zf(i,j)=ber;
        
        % Calculating MMSE filter
        H = fft(h_sym);
%         h_eq = ifft(conj(H)./((abs(H).^2)+(var_n/var_s)));
        h_eq = ifft(conj(H)./((conj(H).*H)));
        data_eq = conv(h_eq,data_r);    
            % Demodulation of MMSE filtered signal
        data_demod_eq = pskdemod(data_eq, M, pi/4);
        [~, ber] = biterr(data_demod_eq(1:L_data), data);
        ber_mmse(i,j)=ber;
        
        %Filtro DFE
        data_eq = dfe(data_r, data_mod);
        
        data_demod_eq = pskdemod(data_eq, M, pi/4);
        [~, ber] = biterr(data_demod_eq(1:L_data), data);
        ber_dfe(i,j)=ber;
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
histogram(ber_mmse(:,7))
title(['Histogram MMSE BER for SNR= ', num2str(SNR(7))],'fontsize', 14);
xlabel('BER', 'fontsize', 12), ylabel('Number of cases', 'fontsize', 12) 