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
SNR = 15;

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


%% Channel adjustment
k = round(rand()*(length(Channel_data.hmat)-1))+1; %Number of CR selected
hmat = Channel_data.hmat;
h_raw = circshift(hmat(:, k), shift); % From all the CR we select a random one
h_raw = h_raw/norm(h_raw); %Normalization of the CR



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

%% Communication

% Generation of the random training data
data = randi([0 1],L_data,1);

% bpsk mapper
data_mod = pskmod(data, M);


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
        H = fft(h_sym);
        h_eq = ifft(1./H);
        data_eq = conv(h_eq,data_r);
    case 2
        disp('Equalization - MMSE');
        disp(['MMSE Equalizer Design: N=', num2str(nTaps)...
            , 'SNR=', num2str(SNR)]);
        % Calculating MMSE filter
        H = fft(h_sym);
        h_eq = ifft(conj(H)./((abs(H).^2)+(var_n/var_s)));
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
    subplot(2,1,1)
    plot(((0:Lsym-1))*q/B*1000,abs(h_sym));
    title('Symbol Spaced Channels Response ','fontsize', 12);
    xlabel('delay [ms]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
    axis([0 40 -inf inf]);
    
    subplot(2,1,2)
    plot(((0:length(h_eq)-1))*q/B*1000,abs(h_eq));
    title('Real','fontsize', 12);
    xlabel('delay [ms]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
    axis([0 40 -inf inf]);
    
    sgtitle('Channel Response','fontsize', 16);
end

scatterplot(data_eq);

data_demod = pskdemod(data_r, M);
data_demod_eq = pskdemod(data_eq, M);

[~, ber] = biterr(data_demod(1:L_data), data)
[~, ber_eq] = biterr(data_demod_eq(1:L_data), data)
