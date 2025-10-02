% AM Modulation and Demodulation 
clc; clear; close all;

%% Parameters
Am = 1;        % Message amplitude
Ac = 2;        % Carrier amplitude
fm = 1000;     % Message frequency (Hz)
fc = 20000;    % Carrier frequency (Hz)
fs = 200000;   % Sampling frequency (Hz)
t = 0:1/fs:2/fm; % Simulation time
mu = Am / Ac;  % Modulation index

%% Message Signal
m = Am * cos(2*pi*fm*t);

%% Carrier Signal
c = Ac * cos(2*pi*fc*t);

%% AM Modulation (DSB-LC)
s = (1 + mu * cos(2*pi*fm*t)) .* c;

%% Add Custom Gaussian Noise
snr = 20; % dB
signal_power = mean(s.^2);
noise_power = signal_power / (10^(snr/10));
noise = sqrt(noise_power) * randn(size(s));
s_noisy = s + noise;

%% Envelope Detection Demodulation
env = abs(s_noisy);

% Simple RC Low-pass Filter implementation
RC = 1/(2*pi*(2*fm)); % cutoff at ~2*fm
dt = 1/fs;
alpha = dt / (RC + dt);
m_rec = zeros(size(env));
for i = 2:length(env)
    m_rec(i) = m_rec(i-1) + alpha*(env(i) - m_rec(i-1));
end

%% Plotting
figure;

subplot(4,1,1);
plot(t, m, 'b');
title('Message Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(4,1,2);
plot(t, c, 'r');
title('Carrier Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(4,1,3);
plot(t, s, 'k');
title('AM Modulated Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(4,1,4);
plot(t, m_rec, 'g');
title('Demodulated Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

%% Spectrum Plot
figure;
nfft = 2^nextpow2(length(s));
f = fs/2*linspace(0,1,nfft/2+1);
S = fft(s, nfft)/length(s);
plot(f, 2*abs(S(1:nfft/2+1)));
title('Spectrum of AM Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
