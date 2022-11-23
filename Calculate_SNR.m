%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate SNR and plot FFT of the Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rSNR = calculate_SNR(SimulationTime,SampleTime,Signal,BW,InSigFreq)

N       = length(Signal);                   % length of the signal vector

%% windowing
f       = 0:floor((N-1)/(2));                 % frequency vector for the plot
f       = f / (N*SampleTime);                         % scale the frequency vector
Signal  = hann(N,'periodic').'.* Signal;                 % hann window
% to lower the leakage effects
%% ToDo: Fast Fourier Transformation
% t = (0:N-1)*SampleTime;
Y = abs(fft(Signal))/N;
%% ToDo: Scaling
%See Windowing
fs = 2*f;
%% ToDo: Limit to F_max = 1/(2*Ts)
signal_fft = 2*abs(Y(1:N/2+1));
figure(1);
plot(fs/fs(length(fs)),pow2db(signal_fft.^2));
xlim([0 0.5]);
title('Input Spectrum (Vin=1/3V; dithering)');
xlabel('frequency/fs');
ylabel('Power(dB)');
%% ToDo: Create window in order to separate signal and noise: around 3-5 Datapoints of input frequency
% SNR depends on Bandwith!
fstep = f(length(f))/(length(f)-1)*2;
if(signal_range >= 2)
    signal_range = [InSigFreq/fstep-2 InSigFreq/fstep+2];
else
    signal_range = [1 InSigFreq/fstep+2];
end
window = ones(1,length(f));
%% ToDo: Separate signal and noise into FFT-Signal and FFT-Noise
signal_window = window * 0;
signal_window(signal_range(1):signal_range(2)) = 1;
signal_separated = signal_fft .* signal_window;
figure(2);
plot(fs, pow2db(signal_separated.^2));
xlim([0 f(length(f))]);
title('Signal Spectrum (Signal @ 100hz)');
xlabel('frequency(Hz)');
ylabel('Power(dB)');

noise_window = window;
noise_window(signal_range(1):signal_range(2)) = 0;
noise_separated = signal_fft .* noise_window;
figure(3);
plot(fs, pow2db(noise_separated.^2));
xlim([0 f(length(f))]);
title('Noise Spectrum (Signal @ 100hz)');
xlabel('frequency(Hz)');
ylabel('Power(dB)');

%% ToDo: Separate noise IN and OUT of bandwidth of 10kH
BW_edge = BW/fstep;

noise_in = noise_separated(1:BW_edge);
figure(4);
plot(fs(1:BW_edge), pow2db(noise_in.^2));
xlim([0 f(length(f))]);
title('In-band Noise Spectrum (Signal @ 100hz)');
xlabel('frequency(Hz)');
ylabel('Power(dB)');
noise_out = noise_separated(BW_edge+1:length(noise_separated));
figure(5);
plot(fs(BW_edge+1:length(noise_separated)), pow2db(noise_out.^2));
xlim([0 f(length(f))]);
title('Out-band Noise Spectrum (Signal @ 100hz)');
xlabel('frequency(Hz)');
ylabel('Power(dB)');
%% ToDo: Calculate SNR
% You can use the snr method in Matlab
% leak            = 50; 
% [~,p]           = max(signal_fft(1:BW_edge));
% sigpos          = [p-leak:p+leak];
% sig_pow         = sum(signal_fft(sigpos).^2); % signal power = sum of magnitudes of bins conrresponding to signal
% signal_fft([sigpos]) = 0; % making all bins corresponding to signal zero:==> what ever that remains is noise 
% noise_pow       = sum(signal_fft.^2); % sum of rest of componenents == noise power
% SNR             = 10*log10(sig_pow/noise_pow);

signal_separated = signal_separated(1:10000);
signal_power = rms(signal_separated);
noise_power = rms(noise_in);
% rSNR = 10*log10(signal_power/noise_power);
rSNR = snr(signal_separated, noise_in);

% signal_power = sum(signal_separated.^2)/5;
% noise_power = sum(noise_in.^2)/length(noise_in);
% rSNR = 10*log10(signal_power/noise_power);
%% ToDo: Plot the spectrum


end