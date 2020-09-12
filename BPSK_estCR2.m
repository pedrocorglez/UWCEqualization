% BPSK UW Acoustic Communication Equalization script
% Author: Pedro Córdoba González
%

close all; clear all;
addpath('Simulated Channel Response'); % We add to the path the folder with the CRs
addpath('Functions'); % We add to the path the folder with the CRs
%% Parameters
%Communication Parameters
M = 2;          % Order of the modulation
L_data = 10000; % Number of transmitted symbols
L_lea = 1000;
Fs_sym = 250;   % Symbol Frequency
SNR = 25;

%Channel Parameters
Fs_h=1e4;       % Sample frequency of Channel Impulse Response
Fs_c=3e4;       % Sample Frequency of Chirp
CRfile='Frequency_Response_sim_dir_45-55kHz_25Hz_60s_0.05s_395_5_25_OK.mat';
Channel_data=load(CRfile); % Data simulated with Stojanovic script
Lf=401; Lt_tot=3603; T_SS=60; T_tot=3*T_SS;
fmin=45e3; % minimum frequency [Hz]
B=10e3; % bandwidth [Hz]
df=25; % frequency resolution [Hz], f_vec=fmin:df:fmax;
dt=50e-3; % time resolution [seconds]
T_SS=60; % coherence time of the small-scale variations [seconds]
shift=10; skip=10;

% Equalizer Parameters
EQ = 1;
nTaps = 20; %Length of the equalizer
delay = 0;
Est = 1;

%% Channel adjustment
k = round(rand()*(length(Channel_data.hmat)-1))+1; %Number of CR selected
hmat = Channel_data.hmat;
H = Channel_data.H_LS;
h_raw = circshift(hmat(:, 83), shift); % From all the CR we select a random one
h_raw = h_raw/norm(h_raw); %Normalization of the CR

% All the CR are shown in 3D plot
figure; axes('fontsize', 16);
% hmat2plot= abs(hmat(1:end, 1:skip:end)).';
s = surf(((0:Lf-1)-shift)/B*1000, (0:skip:Lt_tot-1)*dt, (circshift(abs(hmat(1:end, 1:skip:end)), shift)).',  'CDataMapping','scaled', 'EdgeColor', 'none');
xlabel('delay [ms]', 'fontsize', 16), ylabel('time [s]', 'fontsize', 16), axis('ij');
set(gca,'YTick', 0:T_SS:T_tot);
colorbar;
axis([-1 20 0 180 0 inf]);
view([130.871299093656 39.7200698977097]);
title('Channel Responses through Time');

% The Channel Response is shown
figure;
subplot(3,1,1)
plot(((0:Lf-1)-shift)/B*1000,abs(h_raw));
title('Module','fontsize', 12);
xlabel('delay [ms]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
axis([-inf inf -inf inf]);

subplot(3,1,2)
plot(((0:Lf-1)-shift)/B*1000,real(h_raw));
title('Real','fontsize', 12);
xlabel('delay [ms]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
axis([-inf inf -inf inf]);

subplot(3,1,3)
plot(((0:Lf-1)-shift)/B*1000,imag(h_raw));
title('Imaginary','fontsize', 12);
xlabel('delay [ms]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
axis([-inf inf -inf inf]);

sgtitle('Channel Response','fontsize', 16);

%Resample of the Channel's Response from the first arrival
[p,q] = rat(Fs_sym / Fs_h);
[m,ind] = max(abs(h_raw(1:50)));        %calculation of the first arrival
% h_sym = resample(h_raw(ind:end),p,q);
h_sym = h_raw(ind:q:end);
Lsym = length(h_sym);

q1=q; p1=p;
%% Channel Estimation
% We are going to use a chirp to estimate the channel response as we would
% do in the real communication

%We resample the Channel's Response to fit the chirp
[p,q] = rat(Fs_c / Fs_h);
H_raw = H(:, k);
H_raw_res = resample(H_raw,p,q);
h_raw_res = circshift(ifft(H_raw_res), shift);
% h_raw_res = ifft(H_raw_res);
[m,ind] = max(abs(h_raw_res(1:200)));        %calculation of the first arrival
h_raw_res = h_raw_res(ind:end);
h_raw_res = h_raw_res/norm(h_raw_res); %Normalization of the CR

%we create the chirp to estimate the channel
t = 0:1/Fs_c:1-1/Fs_c; 
swept = chirp(t,0,t(end),B)';

% We conv the two signals
swept_r = cconv(h_raw_res,swept);
swept_r = awgn(swept_r,SNR);

%We calculate the estimated Channel's Response
Y = swept_r(1:Fs_c);
Sxy = conj(swept).*Y;
Sxx = conj(swept).*swept;
Syy = conj(Y).*Y;
H_r = Sxy./Sxx;
h_r = ifft(H_r);

% h_r = circshift(h_r, shift*q/p);
% [m,ind] = max(abs(h_raw_res(1:200)));        %calculation of the first arrival
% h_r_res = h_raw_res(ind:end);
%We resampled to the symbol Frequency
[p,q] = rat(Fs_sym / Fs_c);
h_sym_r = h_r(1:q:end);
h_sym_r=  h_sym_r(1:Lsym);

t = (0:length(h_r)-1)/Fs_c; 
% The Channel Response is shown
figure;
subplot(3,1,1)
plot(t,abs(h_r));
title('Module','fontsize', 12);
xlabel('delay [s]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
axis([-inf 0.04 -inf inf]);

subplot(3,1,2)
plot(t,real(h_r));
title('Real','fontsize', 12);
xlabel('delay [s]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
axis([-inf 0.04 -inf inf]);

subplot(3,1,3)
plot(t,imag(h_r));
title('Imaginary','fontsize', 12);
xlabel('delay [s]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
axis([-inf 0.04 -inf inf]);

sgtitle('Estimated Channel Response','fontsize', 16);

% The Channel Response is shown
figure;
subplot(3,1,1)
plot(t,abs(H_r));
title('Module','fontsize', 12);
xlabel('delay [s]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
% axis([-inf 0.04 -inf inf]);

subplot(3,1,2)
plot(t,real(H_r));
title('Real','fontsize', 12);
xlabel('delay [s]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
% axis([-inf 0.04 -inf inf]);

subplot(3,1,3)
plot(t,imag(H_r));
title('Imaginary','fontsize', 12);
xlabel('delay [s]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
% axis([-inf 0.04 -inf inf]);

sgtitle('Estimated Channel''s Frequency Response','fontsize', 16);

%% Communication

% Generation of the random training data
data = randi([0 1],L_data,1);

% bpsk mapper
data_mod = pskmod(data, M);

% des=randi([0 1],40,1);% desired data symbols
% des=pskmod(des, M);
% rec=conv(h_sym, des); % received data symbols...
% rec = rec(1:40);
% z=filter([1 ;zeros(39,1)],[1],rec);
% Rxx=xcorr(rec);
% Rxz=xcorr(des,z);
% x=toeplitz([Rxx ;zeros(1,1)],zeros(length(Rxx)+1,1));
% h_eq=x\([Rxz; zeros(1,1)]); % coefficients for MMSE equalizer...

% Calculation of the symbols received
data_r_nonoise = conv(h_sym,data_mod);
data_r = awgn(data_r_nonoise,SNR);
var_n = 10^(SNR/10);
var_s = var(data_r_nonoise);

scatterplot(data_r);

switch(EQ)
    case 1
        disp('Equalization - ZFE');
        disp(['ZF Equalizer Design: N=', num2str(nTaps)]);
        % Calculating ZF filter
        H = fft(h_sym_r);
        h_eq = ifft(1./H);
        data_eq = conv(h_eq,data_r);
    case 2
        disp('Equalization - MMSE');
        disp(['MMSE Equalizer Design: N=', num2str(nTaps)...
            , 'SNR=', num2str(SNR)]);
        % Calculating MMSE filter
        H = fft(h_sym_r);
        h_eq = ifft(conj(H)./((abs(H).^2)+(var_n/var_s)));
%         h_eq = ifft(conj(H)./((conj(H).*H)));
        data_eq = conv(h_eq,data_r);
    case 3
        disp('Equalization - DFE');
        dfe = comm.DecisionFeedbackEqualizer('Algorithm','LMS', ...
            'NumForwardTaps',20,'NumFeedbackTaps',10,'StepSize',0.03);
        
        data_eq = dfe(data_r, data_mod);
    otherwise
        disp('EQ no es una opción válida.');
end

if (EQ == 1 | EQ == 2)
    % The Channel Response is shown
    figure;
    subplot(3,1,1)
    plot(((0:Lsym-1))*q/B*1000,abs(h_sym_r));
    title('Absolute','fontsize', 12);
    xlabel('delay [ms]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
    axis([0 40 -inf inf]);
    
    subplot(3,1,2)
    plot(((0:length(h_eq)-1))*q/B*1000,abs(h_eq));
    title('Real','fontsize', 12);
    xlabel('delay [ms]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
    axis([0 40 -inf inf]);
    
    subplot(3,1,3)
    plot(((0:length(h_eq)-1))*q/B*1000,imag(h_eq));
    title('Imaginary','fontsize', 12);
    xlabel('delay [ms]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
    axis([0 40 -inf inf]);
    
    sgtitle('Channel Response','fontsize', 16);
end

scatterplot(data_eq);

data_demod = pskdemod(data_r, M);
data_demod_eq = pskdemod(data_eq, M);

[~, ber] = biterr(data_demod(1:L_data), data)
[~, ber_eq] = biterr(data_demod_eq(1:L_data), data)
