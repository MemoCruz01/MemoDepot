%% Excercise 1
clear; close all; clc;

%% PARAMETERS
Delta = 0.02;              % Window half-width
fs = 50000;                % Sampling frequency (high for accurate FT)
T = 1/fs;
t = -0.05:T:0.05;          % Time axis
N = length(t);

%% SIGNAL DEFINITIONS
% Window w(t)
w = double(abs(t) <= Delta);

% x(t) = 3 sin(100π t + π/4)
x = 3 * sin(100*pi*t + pi/4);

% y(t) = x(t) w(t)
y = x .* w;

%% FOURIER TRANSFORMS (numerical using FFT)
f = (-N/2:N/2-1)*(fs/N);      % frequency axis

W = fftshift(abs(fft(w)));
X = fftshift(abs(fft(x)));
Y = fftshift(abs(fft(y)));

%% ANALYTICAL W(f) FOR RECTANGLE (for reference)
% W(f) = 2Δ sinc( 2Δ f )
W_analytic = 2*Delta * abs(sinc(2*Delta*f));

%% PLOTTING
figure('Name','All Plots','Position',[100 100 1200 900]);

% ----- 1. w(t) -----
subplot(3,2,1);
plot(t, w, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('w(t)');
title('Window w(t)');

% ----- 2. |W(f)| -----
subplot(3,2,2);
plot(f, W_analytic, 'r', 'LineWidth', 2); hold on;
plot(f, W/max(W)*max(W_analytic), 'k--');  % scaled numerical FT
xlim([-300 300]);
xlabel('Frequency (Hz)'); ylabel('|W(f)|');
title('Spectrum of window |W(f)|');

% ----- 3. x(t) -----
subplot(3,2,3);
plot(t, x, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('x(t)');
title('Signal x(t)');

% ----- 4. |X(f)| -----
subplot(3,2,4);
plot(f, X, 'LineWidth', 1.5);
xlim([-300 300]);
xlabel('Frequency (Hz)'); ylabel('|X(f)|');
title('Spectrum |X(f)|');

% ----- 5. y(t) = x(t) w(t) -----
subplot(3,2,5);
plot(t, y, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('y(t)');
title('Windowed signal y(t)=x(t)w(t)');

% ----- 6. |Y(f)| -----
subplot(3,2,6);
plot(f, Y, 'LineWidth', 1.5);
xlim([-300 300]);
xlabel('Frequency (Hz)'); ylabel('|Y(f)|');
title('Spectrum |Y(f)|');

sgtitle('Tutorial: Analog Signals & Spectral Representation (All in One Figure)');




