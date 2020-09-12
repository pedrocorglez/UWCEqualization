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

%Modulation Parameters
M = 4; % Modulation order for QPSK
nfft  = 128;
cplen = 16;
nullIdx = [1:6 65 128-4:128].';
pilotIdx = [12; 26; 40; 54; 74; 88; 102; 116];
numDataCarrs = nfft-length(pilotIdx)-length(nullIdx);
pilots = repmat(pskmod([0 1 2 3 0 1 2 3]',M),1,L_data);
B_mod = 6e3;

%Channel Parameters
Fs_h=1e4;       % Sample frequency of Channel Impulse Response
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
EQ = 3;
nTaps = 20; %Length of the equalizer
delay = 0;


%% Channel adjustment
k = round(rand()*(length(Channel_data.hmat)-1))+1; %Number of CR selected
hmat = Channel_data.hmat;
h_raw = circshift(hmat(:, k), shift); % From all the CR we select a random one
H_raw = Channel_data.H_LS(:,k);
h_raw = h_raw/norm(h_raw); %Normalization of the CR
H_raw = H_raw/norm(H_raw); %Normalization of the CR

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

% The Channel Response is shown
f=fmin:df:fmin+B;
figure;
subplot(3,1,1)
plot(f,abs(H_raw));
title('Module','fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
axis([-inf inf -inf inf]);

subplot(3,1,2)
plot(f,real(H_raw));
title('Real','fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
axis([-inf inf -inf inf]);

subplot(3,1,3)
plot(f,imag(H_raw));
title('Imaginary','fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
axis([-inf inf -inf inf]);

sgtitle('Channel Frequency Response','fontsize', 16);

%We extract the part of the channel which is in the modulation
idc = B/(df*2)+1;   %índice de la frecencia central
id1 = idc - floor(B_mod/(df*2));
id2 = idc + floor(B_mod/(df*2));
H_s = H_raw(id1:id2);

[p,q] = rat(nfft / length(H_s));
H_rs = resample(H_s,p,q);
dataIdx = double(setdiff((1:nfft)',nullIdx));
dataIdx = double(setdiff(dataIdx',pilotIdx));
H = H_rs(dataIdx);

%Resample of the Channel's Response from the first arrival
[p,q] = rat(Fs_sym*nfft / Fs_h);
[m,ind] = max(abs(h_raw(1:50)));        %calculation of the first arrival
h_sym = resample(h_raw(ind:end),p,q);
% h_sym = h_raw(ind:q:end);
Lsym = length(h_sym);

%% Communication

% Generation of the random training data
data = randi([0 M-1],numDataCarrs,L_data);
qpskSym = pskmod(data,M,pi/M);

% dataIn = complex(randn(modDim.DataInputSize));
% pilotIn = complex(rand(modDim.PilotInputSize));
% dataIn = pskmod((randi([0 M-1],modDim.DataInputSize(1),1)),M,pi/M);
% dataIn = randi([0 3],modDim.DataInputSize(1),100);
% pilotIn = ones(modDim.PilotInputSize(1),1);
% pilotIn = randi([0 3],modDim.PilotInputSize(1),100);

% ofdm modulation
data_mod = ofdmmod(qpskSym,nfft,cplen,nullIdx,pilotIdx,pilots);
scatterplot(fft(data_mod))
data_demod = ofdmdemod(data_mod,nfft,cplen,cplen,nullIdx,pilotIdx);
scatterplot(reshape(data_demod,numDataCarrs*L_data,1));

% Calculation of the symbols received
data_r_nonoise = conv(h_sym,data_mod);
data_r = awgn(data_r_nonoise,SNR);
var_s = var(data_r_nonoise);
var_n = 10^((10*log10(var_s)-SNR)/10);


data_r_demod = ofdmdemod(data_r(1:length(data_mod)),nfft,cplen,0,nullIdx,pilotIdx);

scatterplot(reshape(data_r_demod,numDataCarrs*L_data,1));

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
%         h_eq = ifft(conj(H)./((conj(H).*H)+(var_n/var_s)));
        h_eq = ifft(conj(H)./((abs(H).^2)+(var_n/var_s)));
%         H_eq = conj(H)./((abs(H).^2)+(var_n/var_s));
        data_eq = conv(h_eq,data_r);
%         data_eq_demod = repmat(H_eq,1,L_data).*data_r_demod; 
    case 3
        disp('Equalization - DFE');
        dfe = comm.DecisionFeedbackEqualizer('Algorithm','LMS', ...
            'NumForwardTaps',20,'NumFeedbackTaps',10,'StepSize',0.03);
     
        data_eq = dfe(reshape(data_r_demod,[],1), reshape(qpskSym(Llearn),[],1));
        data_eq_demod=reshape(data_eq,size(data_r_demod));
    otherwise
        disp('EQ no es una opción válida.');
end
% 
% if (EQ == 1 | EQ == 2)
%     % The Channel Response is shown
%     figure;
%     subplot(2,1,1)
%     plot(((0:Lsym-1))*q/B*1000,abs(h_sym));
%     title('Symbol Spaced Channels Response ','fontsize', 12);
%     xlabel('delay [ms]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
%     axis([0 40 -inf inf]);
%     
%     subplot(2,1,2)
%     plot(((0:length(h_eq)-1))*q/B*1000,abs(h_eq));
%     title('Real','fontsize', 12);
%     xlabel('delay [ms]', 'fontsize', 12), ylabel('Amplitude', 'fontsize', 12)
%     axis([0 40 -inf inf]);
%     
%     sgtitle('Channel Response','fontsize', 16);
% end
if EQ ~= 3
[data_eq_demod,pilots] = ofdmdemod(data_eq(1:length(data_mod)),nfft,cplen,cplen,nullIdx,pilotIdx);
end

scatterplot(reshape(data_eq_demod,[],1));

data_demod = pskdemod(data_r_demod, M);
data_demod_eq = pskdemod(data_eq_demod, M);

[~, ber] = biterr(data_demod, data)
[~, ber_eq] = biterr(data_demod_eq, data)