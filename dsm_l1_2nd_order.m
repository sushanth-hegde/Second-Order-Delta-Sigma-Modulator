%% Init
clear all;
close all;

%% Path
addpath('../delsig')

%% Spec.
OSR = 64;     % simulation length (output samples)
order= 2;     % number of sinusoids
N = 2^16;
fN = 44.1e3;  % Nyquist freq
fs = fN*OSR;  % sampling frequency 
Ts  = 1/fs;   % time step
fbin = fs/N;
nlev = 2;
%% Input tone
tone = 10e3;                % 10k
A  = 0.8;                   % signal amplitude fi 
offset = 0;
fx = ceil(tone/fbin) * fbin;% input tone fbin
t=0:Ts:N*Ts;
u= A*sin(2*pi*fx*t);

%% Coefficients
form = 'CIFB';
ntf = synthesizeNTF(order,OSR); 
[a,g,b,c] = realizeNTF(ntf,form);

% poles plot NTF,STF plot
plotPZ(ntf,'b',10,0);
b(2:end) = 0; % for a maximally flat STF
ABCD = stuffABCD(a,g,b,c,form);
[Ha Ga] =  calculateTF(ABCD);%NTF and STF
f = linspace(0,0.5,100);
z = exp(2i*pi*f);
magHa  = dbv(evalTF(Ha,z));
magGa = dbv(evalTF(Ga,z));
figure,plot(f,magHa, 'b', f,magGa,'--');
title('2nd order')
title('NTF and STF plot for the modulator realization')
ylabel('Magnitude in dB')
xlabel('normalized frequency')
vsim = simulateDSM(u,ntf);
%% Simulation
mdl = 'dsm_l1_sim';
load_system(mdl);
open_system(mdl);

simoptions = simset('Solver','FixedStepDiscrete', ...
    'RelTol',1e-3, ...
    'MaxStep', Ts);

[t_, u_, vsim] = sim(mdl, max(t), simoptions, [t',u']);

Q = vsim(:,1);
Y = vsim(:,2);
 %% Psd

N = length(vsim(:,2));
dft = fft(vsim(:,2));
dft = dft(1:N/2+1);
psdx = (1/(fs*N)) * abs(dft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
yy = smooth(psdx)
freq = 0:fs/length(vsim(:,2)):fs/2;

figure;
plot(freq,10*log10(yy))
grid on
title ('FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequncy (dB/Hz)')

%% Post-processing
figure(1)
plot(t, u, 'rx-', t, Q, 'g');
axis([ min(t) max(t) 1.1*min(Q) 1.1*max(Q) ]);
hold on;
stairs(t, Y, 'b');
xlabel('Time   [ t/T ]');
ylabel('Amplitude');
legend('X', 'Q', 'Y');
title('2nd Order Sigma-Delta');
hold off;
%%  peak snr

[snr_pred,amp] = predictSNR(ntf,OSR);
[snr,amp] = simulateSNR(ntf,OSR);
plot(amp,snr_pred,'b',amp,snr,'gs');
grid on;
figureMagic([-100 0], 10, 2,[0 100], 10, 1);
xlabel('Input Level, dB');
ylabel('SNR dB');
s=sprintf('Peak SNR = %4.1f dB\n',max(snr));
s1=sprintf('Predicted SNR = %4.1f dB\n',max(snr_pred));
Rbit = ((max(snr))-1.76)/6.02; % Equivalent resolution oin bits
text(-35,15,s);
text(-35,20,s1);
text(-35,25,sprintf('ENOB = %4.1f Bits\n',Rbit));

%% Decimation
%%sinc3 filter
N=8;
Fs=2.822e06;
h1 = ones(1,N)./N;
h2 = conv(h1,h1);
h3 = conv(h1,h2);

sinc = dsp.FIRDecimator(N,h3);
y = sinc(vsim(1:end-1,2));
% fvt = fvtool(y,'Fs',Fs,'Color','white');
% fvt.MagnitudeDisplay = 'Zero-phase';
% fvt.Analysis = 'impulse';
%% halfband

Fs  = 352.75e3;
Fp  = 20e3;

num = firhalfband(10,Fp/(Fs/2));
% fvt = fvtool(num,'Fs',Fs,'Color','white');
% fvt.MagnitudeDisplay = 'Zero-phase';
% fvt.Analysis = 'impulse';

 FrameSize = 256;
FsIn = Fs;
 N   = 6;
halfbandDecimator = dsp.FIRHalfbandDecimator('SampleRate',FsIn, ...
    'Specification','Filter order and transition width', ...
    'FilterOrder',N,'TransitionWidth',100);
% fvtool(halfbandDecimator,'Fs',FsIn,'Color','white');
% scope = dsp.SpectrumAnalyzer('SampleRate',Fs,'SpectralAverages',5);
    y1 = halfbandDecimator(y);      %176KHz                
    xd = halfbandDecimator(y1);     %88KHz
    xd1 = halfbandDecimator(xd);    %44.1KHz(nyquist rate)
   %scope(xd);
%release(scope);



%% fft
N = length(xd1);
dft = fft(xd1);
dft = dft(1:N/2+1);
psdx = (1/(fN*N)) * abs(dft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
yy = smooth(psdx)
freq = 0:fN/length(xd1):fN/2;
figure;
plot(freq,10*log10(psdx))
grid on
title ('FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequncy (dB/Hz)')
