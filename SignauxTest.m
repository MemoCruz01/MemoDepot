%% Exam 1 Signaux Part 1
clear; close all; clc;


% 1.-Let w(t) be 1 if t is from a to b; 0 otherwhise.
% Determine W(f)=FourrierTransform{w(t)}.
% 
% 2.-Let x(t) = cos(t) Determine X(f)=FourrierTransform{x(t)}.
% 3.-Let y(t) = exp(-t) Determine Y(f)=FourrierTransform{y(t)}.
% 
% 4.- Find the expression of the signal h(t)= x(t) w(t) y(t)
% Then compute Fourrier Transform
% 
% Write a Matlab code doing the 4 steps and plot on a same graph the spectra signals x(t) and h(t). Conclude
% on your observations

% PARAMETERS
a = 0; 
b = 5;                  % time window of w(t)
Fs = 2000;              % sampling frequency
T = 1/Fs;               
t = a:T:b;              % time vector
N = length(t);

% 1) Rectangular window w(t)
w = ones(size(t));      % from a to b

% 2) x(t) = cos(t)
x = cos(t);

% 3) y(t) = exp(-t)
y = exp(-t);

% 4) h(t) = x(t) w(t) y(t)
h = x .* w .* y;

% FREQUENCY AXIS
f = (-N/2:N/2-1)*(Fs/N);

% Fourier transforms (numerical)
W = fftshift(fft(w));
X = fftshift(fft(x));
Y = fftshift(fft(y));
H = fftshift(fft(h));

% NORMALIZE for visualization
W = abs(W)/max(abs(W));
X = abs(X)/max(abs(X));
Y = abs(Y)/max(abs(Y));
H = abs(H)/max(abs(H));

% PLOTS
figure;
plot(f, X, 'LineWidth', 1.5); hold on;
plot(f, H, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectra of x(t) and h(t)');
legend('X(f)', 'H(f)');
grid on;

% Conclusion
% Conclusion / Interpretation
% 
% The spectrum X(f) consists of two impulses at
% 
% f= ±2π1 ​≈ ±0.159 Hz.
% 
% The signal h(t) multiplies x(t) by 2 proportions:

% -an exponential decay, broadening the spectrum
% -a rectangular window → creates sinc-like spreading, introducing many new frequency components.




%% Exam 1 Signaux Part 2

% 1.- Condition to be fulfilled is the Nyquist-Shannon. When sampling a signal into the discrete domain
% in order to be able to reconstruct / retrieve information of the main signal,
% the sampling frequency must be at least as twice of the maximum freauency of the main signal.
%             Fs >= 2 * Fmax                      In practical scenarios is usually Fs > 2 * Fmax
% 
% for x(t)=cos(t)=cos(2pi*1/(2pi)t) the max frequency is f=1/(2pi) Hz
% so the sampling frequency must be at least Fs=2 * 1/(2pi) =1/(pi) Hz to ful fill Shannon Theorem
% 
% 
% 2.- Condition for h(t) to be satisfied:    what would be the sampling frequency to convert h(t) to a digital signal h_n. Shannon Theory
% 
% For this explession it will still be Fs=1/(pi), as the function y(t) = e^(-t) is a low pass filter
% which is not providing nez higher frequencies.

% 3.- Provide relationship between h(t) with h_n.
% 
% The relqtionships rely on x[n]=x(nT) and h[n]=h(nT). Being Ts the sampling period and Fs the Sampling Freq;
% Fs=1/Ts
% 
% Therefore using x[n]=x(nT), we get in discrete part:
% x[n]=cos(nTs​)
% 
% Therefore using h[n]=h(nT), we get in discrete part:
% h[n]=cos(nTs​) * e^(-nTs) * ​w(nTs​)

% 4.- For a given frequency Fs=10Hz, whats the number of samples N required to correctly represent x(t).
% 
% Fundamental period of x(t)          x(t): T0=2pi
% Numerically:   T0=6.283s
% ​
% To represent one full period:
% N=Fs​*T0​ = 10 * 6.283 = 62.83
% 
% Aprox 63 samples are required for 10Hz
	​
%% 5.- Matlab: Consider N=8. Determine the Discrete Fourrier Transform of x_n and h_n. Plot the spectrum of these signals and conclude
clear; close all; clc;
% PARAMETERS
N = 8;              % number of samples
Fs = 10;            % sampling frequency
Ts = 1/Fs;
n = 0:N-1;          % discrete time indices

% SIGNALS
x_n = cos(n*Ts);               % sampled cosine
h_n = x_n .* exp(-n*Ts);       % simple version without window for clarity

% DFT
X_k = fft(x_n);
H_k = fft(h_n);

% Frequency axis
k = 0:N-1;
f = k*Fs/N;

% PLOTS
figure;
subplot(2,1,1);
stem(f, abs(X_k), 'LineWidth', 1.5);
title('|X[k]| - DFT of x[n]');
xlabel('Frequency (Hz)');
ylabel('Magnitude'); grid on;

subplot(2,1,2);
stem(f, abs(H_k), 'LineWidth', 1.5);
title('|H[k]| - DFT of h[n]');
xlabel('Frequency (Hz)');
ylabel('Magnitude'); grid on;





