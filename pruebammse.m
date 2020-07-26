clc
close all
clear all
H = comm.QPSKModulator('BitInput',true);
Hdemod = comm.QPSKDemodulator('BitOutput',true);
% hScope = commscope.ScatterPlot;
% hScope.Constellation = [0.7071+0.7071i -0.7071+0.7071i -0.7071-0.7071i 0.7071-0.7071i];
% hScope.SamplesPerSymbol = 1;
n=16;    %no of random data...
r_data = randi(n,1); %random numbers generator...
data_qpsk=[];
data_qpsk = step(H,r_data); %Data converted in QPSK symbols...
d1=data_qpsk(1);
d5=data_qpsk(5);
data_qpsk(1) = 0.7071+0.7071i; % adding two pilots at location 1 and 5.
data_qpsk(5) = 0.7071+0.7071i;
% update(hScope, data_qpsk);
dawgn=awgn(data_qpsk,0); % Adding white Gaussian Noise
est(1)=dawgn(1);
est(2)=dawgn(5);
% MMSE starts here......
des=[0.7071+0.7071i 0.7071+0.7071i];% desired data symbols
rec=[est(1) est(2)]; % received data symbols...
z=filter([1 0 0 0 0 0 0 0],[1],rec);
Rxx=xcorr(rec);
Rxz=xcorr(des,z);
x=toeplitz([Rxx zeros(1,5)],zeros(1,8))
cof=x\([Rxz zeros(1,5)].'); % coefficients for MMSE equalizer...
det1=filter(cof,[1],dawgn);
det1
for i=1:8
det(i)=filter(cof,[1],dawgn(i));
end
det.'
