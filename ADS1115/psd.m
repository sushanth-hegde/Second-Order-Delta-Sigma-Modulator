close all;
clear all;
Fs=51;
X = xlsread('b1.xlsx');
x = X(:,1);
nfft = 2^nextpow2(length(x));
Pxx = abs(fft(X,nfft)).^2/length(x)/Fs;
%to Create a single-sided spectrum
Hpsd = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);
plot(Hpsd);


