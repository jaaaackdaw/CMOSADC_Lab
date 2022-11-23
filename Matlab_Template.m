%% --------------------------MODULATOR SYNTHESIS---------------------------
% 
% This exemplary code allows to quickly design CT modulators by utilizing the 
% functions of the Schreier's Delta-Sigma Toolbox. It works only for Low-Pass
% modulators and Butterworth NTFs approximations. Other types of
% Sigma-Delta modulators would require modifications of the code.
%
% To make your own design simply change the parameters in STEP[1] only,
% unless you can confidently code in Matlab and know the Schreier's Toolbox
% functions available.
%
% Note that this code does not map with any of the SIMULINK models of the
% Toolbox.
%
% ALL functions used in this code are covered by:
% Copyright (c) 2009, Richard Schreier - All rights reserved.

clc
close all

%% STEP [1] - MODULATOR SPECIFICATIONS:
% ToDo: Rewrite all NaNs!
order = 2; %1%2     % Order of the modulator
OSR = 256; %21      % Oversampling ratio of the modulator 
form = 'FB';    % Modulator architecture (accepted entries include: 
                    % 'FF' and 'FB' for feedforward and feedback structures 
                    % respectively)
form1 = 'CIFB';
nLev = 2;       % Number of levels in the quantizer 
OBG = 1.5;      % Out of band gain of NTF - NaN for default by Schreier's
                    % toolbox(1.5 - Lee criterion)
opt = 1;        % Optimization of NTF zeroes (0 = no opt, 1 = opt) 
f0 = 0;         % Centre frequency - it is advised to leave zero for CT Low
                    % Pass modulators unlessfully confindent with CT-MOD behavior                    
tdac = [0 1];	% DAC timing. [0 1] means zero-delay non-return-to-zero
                    % other timing may include [0 0.5] for RZ pulses, [0.5 1]
                    % for HRZ, etc.
N=NaN;         % N-point FFT
BW = 1e4;       % Band width

frequency = 0;
%Simulink variables - you can use these variables directly in simulink
SampleFreq = OSR*2*BW;
SampleTime = 1/SampleFreq;
SimulationTime = 1;
InSigAmp = 1/3;
InSigFreq = correctFrequToSim(frequency, SimulationTime);
                   
%% ToDO: STEP [2] - SYNTHESIZE NTF:

if exist('LiveDemo','var') == 0
    LiveDemo=0;
end
format compact;
fig1pos1 = [9 630 200 200];
fig1pos2 = [10 407 480 420];
% fig2pos1 = [239 630 450 200];
% fig2pos2 = [241 341 523 485];

clc
fprintf(1,'\t\tNTF Synthesis-- 2nd-order modulator\n\n');
echo on
% order = 2;
% OSR = 128;
% opt = 1;
H = synthesizeNTF(order,OSR,opt);
H
echo off

figure(1); clf
plotPZ(H);
set(gcf,'MenuBar','none'); 
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','NTF Poles and Zeros');
if LiveDemo
    set(1,'position',fig1pos2);
    changeFig(18,2,12);
else
	fprintf(1,'paused\n');
end
pause
if LiveDemo
    set(1,'position',fig1pos1);
    changeFig;nLev
end

figure(2); clf
f = [linspace(0,0.75/OSR,100) linspace(0.75/OSR,0.5,100)];
z = exp(2i*pi*f);
magH = dbv(evalTF(H,z));
subplot(211);
plot(f,magH);
figureMagic([0 0.5],0.05,2, [-100 10],10,2 );
xlabel('Normalized frequency (1\rightarrow f_s)');
ylabel('dB');
set(gcf,'MenuBar','none'); 
set(gcf,'NumberTitle','off'); 
set(gcf,'Name','NTF Magnitude Response');
fstart = 0.01;
f = linspace(fstart,1.2,200)/(2*OSR); z = exp(2i*pi*f);
magH = dbv(evalTF(H,z));
subplot(212);
semilogx(f*2*OSR,magH);
axis([fstart 1.2 -100 -20]);
grid on
sigma_H = dbv(rmsGain(H,0,0.5/OSR));
hold on;
semilogx([fstart 1], sigma_H*[1 1]);
plot([fstart 1], sigma_H*[1 1],'o');
text( 0.15, sigma_H+5, sprintf('rms gain = %5.0fdB',sigma_H));
xlabel('Normalized frequency (1\rightarrow f_B)');
ylabel('dB')
if LiveDemo
    set(2,'position',fig2pos2);
    changeFig(18,2,12);
else
	fprintf(1,'paused\n');
end
% pause
if LiveDemo
    set(2,'position',fig2pos1);
    changeFig;
    
end

%% Simulation DSM
% order = 2;
% OSR = 32;
% opt = 1;
% nLev = 2;
% H = synthesizeNTF(order,OSR,opt);

% Nfft = 2^13;
% tone_bin = 100;
t=(0:1/(OSR*2*BW):1);
u=0.5*sin(2*pi*BW*t);
v=simulateDSM(u,H,nLev);
n=1:OSR*2;
stairs(t(n),v(n),' b');
hold on;
stairs(t(n),u(n),' g');
title('Input Output Simulation of DSM');
xlabel('Time'); ylabel('Amplitude');

%% ToDO: STEP [3] - CONTINUOUS TIME MAPPING (State-Space method):
% Realize the CT NTF from 
% Partition the coefficient for ABCD matrix to achieve coefficient matrix
[a,g,b,c] = realizeNTF(H, form1)
ABCD = stuffABCD(a,g,b,c,form1)
[ABCDs,umax] = scaleABCD(ABCD,nLev,f0,0.9)
[a,g,b,c] = mapABCD_CIFB_2nd_LSE(ABCDs)
% [a,g,b,c] = mapABCD(ABCDs,form1)
% Order=2;OSR=32;nLev =2;opt=1;
% Xlim=0.9;f=0;H = synthesizeNTF (order,OSR,opt);
% Form = 'CIFB';
% [a,g,b,c] = realizeNTF (H,Form)
% ABCD = stuffABCD (a,g,b,c,Form);
% [ABCDs,umax]=scaleABCD (ABCD,nLev,f,Xlim); 
% [a,g,b,c] = mapABCD_CIFB_2nd_LSE(ABCDs)
% [a,g,b,c] = mapABCD(ABCDs,Form) %scaled ABCD matrix

%% STEP [4] - DYNAMIC RANGE SCALING 
% ToDo: Use simulateDSM to find xmax

% ToDo: Scaling with singal limits
% Scale your matrix with Schreier

% ToDo: Map coefficients with LSE methode to be consistent with task
% description

%% Simulink plot
% Simulate your simulink model with the calculated coefficients and save
% the signal output
out = sim('second_order_dithering');
h = (out.simout.Data).';
% plot spectrum and calculate SNR of the output signal with your
% calculate_SNR - methode
SNR = Calculate_SNR(SimulationTime,SampleTime,h,BW,InSigFreq);