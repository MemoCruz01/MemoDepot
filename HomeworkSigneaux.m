%% HomeworkSigneaux 

clear; close all; clc;

%Professor, I ,odified the code to use different parameters. First of all
%depending on the sampling frequency; the aliasing will exist or not, so
% modifying my sample freq / 2 greater than 44.1kHz allowed me to she thepresence of the aliasing.
% Then, modifying the values of A,B, a, b; c, I saw chenges on the
% amplitude and the frequency at which the signal shifts. Cause the sin
% func depends on time, a, b and c will modify this rate. Then A and B
% directly modify the amplitude. Also A when 0 allows to see the behavior
% of the pure cosine (chirp disappears produced by sine). And ACF goes
% periodically

%% -------------------------------------------------------
% PARAMETERS (from your table)
A = 1;
B = 2;
a = 1;
b = 0.5;
c = 1;
f0 = 5000;                % Hz

%% -------------------------------------------------------
% TIME AXIS AND SAMPLING CONDITION
fs = 90555; %44100;               % Standard audio sampling rate
T = 3;                    % Duration (seconds)
t = 0:1/fs:T-1/fs;

% Define signal
x = (A/pi)*sin(2*pi*f0*(a*t.^2 + b*t + c)) + ...
    (B/(3*pi))*cos(2*pi*f0*t);

%% -------------------------------------------------------
% 1. SAMPLING CONDITION
fprintf("\n---- 1. SAMPLING CONDITION ----\n");
max_inst_freq = f0 * (2*a*t(end) + b);  % instantaneous frequency: f(t)=f0*(2at+b)
fprintf("Maximum instantaneous frequency = %.2f Hz\n", max_inst_freq);
fprintf("Nyquist frequency = fs/2 = %.2f Hz\n", fs/2);

if max_inst_freq < fs/2
    fprintf("=> Sampling IS sufficient (no aliasing expected)\n");
else
    fprintf("=> Sampling is NOT sufficient (aliasing will occur!)\n");
end

%% -------------------------------------------------------
% 2. SPECTRUM OF x(t)
N = length(x);
X = fftshift(abs(fft(x)));
f = (-N/2:N/2-1)*(fs/N);

%% -------------------------------------------------------
% 3. ADD NOISE
x_noisy = x + 0.2*randn(size(x));  % zero-mean Gaussian noise
X_noisy = fftshift(abs(fft(x_noisy)));

%% -------------------------------------------------------
% 4. MEAN VALUE (global + local means at 3 instants)

% Global mean
mean_x = mean(x);
fprintf("\n---- 4. MEAN VALUE ----\n");
fprintf("Global mean of x(t) = %.6f\n\n", mean_x);

% Ask user for 3 time instants
%fprintf("Enter THREE time instants inside [0, %.3f] seconds:\n", t(end));
t1 = 1; %input('t1 = ');
t2 = 2; %input('t2 = ');
t3 = 2.5; %input('t3 = ');

% Local mean parameters
window_length = 0.05;    % 50 ms window (adjustable)
Nw = round(window_length * fs);

% Function to compute local mean
local_mean = @(tc) mean( x(abs(t - tc) <= window_length/2) );

% Compute local means
m1 = local_mean(t1);
m2 = local_mean(t2);
m3 = local_mean(t3);

% Display results
fprintf("\nLocal Mean Values (using %.0f ms window):\n", window_length*1000);
fprintf("Mean around t1 = %.4f → %.6f\n", t1, m1);
fprintf("Mean around t2 = %.4f → %.6f\n", t2, m2);
fprintf("Mean around t3 = %.4f → %.6f\n", t3, m3);


%% -------------------------------------------------------
% 5. AUTOCORRELATION
[Rx, tau] = xcorr(x, 'normalized');
[Rx_noisy, ~] = xcorr(x_noisy, 'normalized');

%If A = 0, the chirp disappears and only cosine remains → autocorr becomes periodic.

%% -------------------------------------------------------
% 6. EXPLORATION: spectrogram
window = hamming(1024);
noverlap = 512;
nfft = 2048;


% %% Some explanations
% Because this frequency depends on time (there is 𝑡 t inside), the signal is a chirp.
% A spectrogram shows:
% 
% Time on X-axis
% Frequency on Y-axis (you enabled 'yaxis')
% Color/brightness = amplitude at that time and frequency
% So, when the signal frequency increases over time, the spectrogram shows a slanted bright line going upward.
% This straight upward line = linear chirp.
%% -------------------------------------------------------
% PLOTS — ALL IN ONE FIGURE
figure('Name','Lab Session Results','Position',[100 50 1200 900]);

% ---- 1. x(t)
subplot(3,2,1);
plot(t(1:2000), x(1:2000));  % short window to show oscillations
title('Original Signal x(t)');
xlabel('Time (s)'); ylabel('Amplitude');

% ---- 2. Spectrum
subplot(3,2,2);
plot(f, X); xlim([-20000 20000]);
title('Spectrum of x(t)');
xlabel('Frequency (Hz)'); ylabel('|X(f)|');

% ---- 3. Noisy signal
subplot(3,2,3);
plot(t(1:2000), x_noisy(1:2000));
title('Noisy Signal x_{noisy}(t)');
xlabel('Time (s)'); ylabel('Amplitude');

% ---- 4. Spectrum with noise
subplot(3,2,4);
plot(f, X_noisy); xlim([-20000 20000]);
title('Spectrum of Noisy Signal');
xlabel('Frequency (Hz)'); ylabel('|X(f)|');

% ---- 5. Autocorrelation
subplot(3,2,5);
plot(tau/fs, Rx); hold on;
plot(tau/fs, Rx_noisy, 'r');
legend('Clean','Noisy');
title('Autocorrelation');
xlabel('Lag (s)'); ylabel('R_x(\tau)');

% ---- 6. Spectrogram
subplot(3,2,6);
spectrogram(x, window, noverlap, nfft, fs, 'yaxis');
title('Spectrogram of x(t)');

sgtitle('Complete Lab Session Signals and Analysis');


