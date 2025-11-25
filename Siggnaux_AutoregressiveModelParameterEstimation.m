%% ================================================================
%  AR(5) PARAMETER ESTIMATION USING AUTOCORRELATION (YULE–WALKER)
%  Dataset: ar5_dataset.csv
%  Author: Gus (polished by ChatGPT)
% ================================================================

clear; clc; close all;

%% --- 1. Load dataset ---
% The CSV must contain a single column of time-series values.
data = readmatrix('C:\Users\o22508808\Downloads\ar5_dataset.csv');
x = data(:);     % ensure column vector
N = length(x);

fprintf('Loaded %d samples.\n', N);

%% --- 2. Choose AR model order ---
p = 5; %3 to 7
%Number of parameters

%% --- 3. Autocorrelation estimation ---
% Use unbiased estimation for correct YW formulation.
[r_full, lags] = xcorr(x, p, 'unbiased');
% Keep r(0) ... r(p)
r = r_full(lags >= 0);

%% --- 4. Build Toeplitz autocorrelation matrix R ---
R = toeplitz(r(1:p));      % R(i,j) = r(|i-j|)
r_vec = r(2:p+1);          % [r1 r2 ... rp]^T

%% --- 5. Solve Yule–Walker: R*a = -r_vec ---
a = -R \ r_vec            % AR parameters [a1..ap]

%% --- 6. Estimate white-noise variance ---
sigma2_e = r(1) + a' * r_vec;

%% --- 7. Compare with MATLAB built-in aryule() ---
[a_yw, sigma2_yw] = aryule(x, p)

%% --- 8. Compute residual error of Yule–Walker model ---
% e[k] = x[k] + a1 x[k-1] + ... + ap x[k-p]
e = filter([1 a'], 1, x);

%% --- 9. Display results ---
fprintf('\n===== AUTOCORRELATION-BASED (YULE–WALKER) RESULTS =====\n');
fprintf('Estimated AR coefficients (manual):\n');
disp(a');
fprintf('Noise variance estimate:  %.6f\n', sigma2_e);

fprintf('\nMATLAB aryule() coefficients:\n');
disp(a_yw(2:end));
fprintf('MATLAB aryule variance:   %.6f\n', sigma2_yw);

%% --- 10. Plot residuals ---
figure;
subplot(2,1,1);
plot(e); grid on;
title('Residual Error (Yule–Walker)');
xlabel('Sample'); ylabel('e[k]');

subplot(2,1,2);
autocorr(e);
title('Autocorrelation of Residuals');

%% ======================================================
%  Plot the raw AR(5) dataset
%  File: ar5_dataset.csv
% ======================================================
% 
% % --- Plot the signal ---
% figure;
% plot(x, 'LineWidth', 1.2);
% grid on;
% 
% title('Raw Time-Series Data from ar5\_dataset.csv');
% xlabel('Sample Index k');
% ylabel('x[k]');
% 
% % Optional: zoom into first 200 samples if the dataset is long
% figure;
% plot(x(1:200), 'LineWidth', 1.2);
% grid on;
% title('First 200 Samples (Zoomed View)');
% xlabel('Sample Index k');
% ylabel('x[k]');



%% ======================================================
% 11. Produce AR model prediction and compare with data
% ======================================================

% Preallocate predicted signal
x_hat = zeros(size(x));

% Start from k = p+1 because model needs p past values
for k = p+1:N
    x_hat(k) = -a' * x(k-1:-1:k-p);
end

% Plot comparison
figure;
plot(x, 'LineWidth', 1.2); hold on;
plot(x_hat, 'r', 'LineWidth', 1.2);
grid on;

legend('Original data x[k]', 'AR model prediction \hat{x}[k]');
title(sprintf('AR(%d) Model — Prediction vs Raw Data', p));
xlabel('Sample index k');
ylabel('Amplitude');



%% ======================================================
% 13. Analyze residual distribution using dfittool
% ======================================================

% Use only valid residuals
e_valid = e(p+1:end);

% Launch distribution fitting tool
dfittool(e_valid);

% Optional: print basic stats
fprintf('\n=== Residual Basic Statistics ===\n');
fprintf('Mean:     %.6f\n', mean(e_valid));
fprintf('Variance: %.6f\n', var(e_valid));
fprintf('Skewness: %.6f\n', skewness(e_valid));
fprintf('Kurtosis: %.6f\n', kurtosis(e_valid));
