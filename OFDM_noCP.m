% BPSK UW Acoustic Communication Equalization script
% Author: Pedro Córdoba González
%

close all; clear all;
addpath('Simulated Channel Response'); % We add to the path the folder with the CRs
addpath('Functions'); % We add to the path the folder with the CRs
%% Parameters
%Communication Parameters
L_data = 1; % Number of transmitted symbols
L_sym= 128;
L_lea = 1000;
Fs_sym = 250;   % Symbol Frequency
SNR = 15;

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

figure(1)
scatter(dataCarriers,zeros(1,length(dataCarriers)),'filled')
hold on
scatter(pilotCarriers,zeros(1,length(pilotCarriers)),'filled')
axis([0 K -1 1])
grid on
legend('data', 'pilots')

%% Communication

% We create the bits for each symbol
dataL_sym = length(dataCarriers); % number of payload bits per OFDM symbol
data = randi([0 M-1],dataL_sym,1); %K random data 

% We modulate the data in PSK symbols
data_psk = pskmod(data, M, phase);
pilot_psk = pskmod(zeros(P,1),M, phase); % The known value each pilot transmits

symbol(dataCarriers) = data_psk;
symbol(pilotCarriers) = pilot_psk;

figure(2)
stem(dataCarriers, real(symbol(dataCarriers)),'b')
hold on
grid on
stem(pilotCarriers, real(symbol(pilotCarriers)),'r')

% We calculate the symbol in time 
symbol_t = ifft(symbol);

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

% Channel convolution
symbol_t_r = conv(symbol_t, h);
symbol_t_r = symbol_t_r(1:K);
symbol_t_r = awgn(symbol_t_r, SNR, 'measured');

% Back to Frequency domain
symbol_r = fft(symbol_t_r);

% Channel estimation
pilots_r = symbol_r(pilotCarriers); %Extraction of the pilots from the received symbol
H_est_pilots = pilots_r ./ pilot_psk; % divide by the transmitted pilot values

x = 1:P;
xq = 1:(P-1)/K:P-(P-1)/K;
%interpolación lineal
H_est_a = interp1(x, abs(H_est_pilots), xq, 'linear');
H_est_p = interp1(x, angle(H_est_pilots), xq, 'linear');
H_est = H_est_a .* exp(1i*H_est_p);
%interpolación cuadrática
H_est_a = interp1(x, abs(H_est_pilots), xq, 'spline');
H_est_p = interp1(x, angle(H_est_pilots), xq, 'spline');
H_est_2 = H_est_a .* exp(1i*H_est_p);

figure(4)
plot(allCarriers, real(H),'g')
hold on
grid on
stem(pilotCarriers, real(H_est_pilots))
plot(allCarriers, real(H_est),'b')
plot(allCarriers, real(H_est_2),'r')
legend('Correct','Pilots','Linear int.','Cuadratic int.')
xlim([0 K+1])

figure(5)
plot(allCarriers, abs(H),'g')
hold on
grid on
stem(pilotCarriers, abs(H_est_pilots))
plot(allCarriers, abs(H_est),'b')
plot(allCarriers, abs(H_est_2),'r')
legend('Correct','Pilots','Linear int.','Cuadratic int.')
xlim([0 K+1])

% Equalization
% symbol_eq = symbol_r ./ H_est';
% symbol_eq_2 = symbol_r ./ H_est';

% symbol_eq = (conj(H_est').*symbol_r) ./ abs(H_est').^2;
% symbol_eq_2 = ( conj(H_est').*symbol_r )./ abs(H_est').^2;

symbol_eq = (conj(H).*symbol_r) ./ abs(H).^2;
symbol_eq_2 = ( conj(H).*symbol_r )./ abs(H).^2;

% Demodulation without equalization
data_r = pskdemod(symbol_r(dataCarriers), M, phase);
[~, ber] = biterr(data_r, data)

% Demodulation with equalization
data_r_eq = pskdemod(symbol_eq(dataCarriers), M, phase);
[~, ber_eq] = biterr(data_r_eq, data)

% Demodulation with equalization
data_r_eq2 = pskdemod(symbol_eq_2(dataCarriers), M, phase);
[~, ber_eq] = biterr(data_r_eq2, data)

symbol_eq_n = (symbol_eq(dataCarriers)-mean(symbol_eq(dataCarriers)))/std(symbol_eq(dataCarriers))*0.7;
symbol_r_n = (symbol_r(dataCarriers)-mean(symbol_r(dataCarriers)))/std(symbol_r(dataCarriers))*0.7;


% We show the results
figure
subplot(1,2,1)
X=[real(data_psk) real(symbol_r_n)];
Y=[imag(data_psk) imag(symbol_r_n)];
plot(X', Y','b')
hold on
grid on
scatter(X(:,1),Y(:,1),'filled','r')
scatter(X(:,2),Y(:,2),'filled','b')
legend('deviation', 'data', 'data rec')
title('Symbols received');

subplot(1,2,2)
X=[real(data_psk) real(symbol_eq_n)];
Y=[imag(data_psk) imag(symbol_eq_n)];
plot(X', Y','b')
hold on
grid on
scatter(X(:,1),Y(:,1),'filled','r')
scatter(X(:,2),Y(:,2),'filled','b')
legend('deviation', 'data', 'data rec')
title('Symbols received after equalization');

