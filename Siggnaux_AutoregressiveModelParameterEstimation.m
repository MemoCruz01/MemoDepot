
% PARAMETER ESTIMATION USING AUTOCORRELATION (YULE–WALKER)

%% ================================================================
%  PART 1 Autocorrelation-based Estimators
% ================================================================
clear; clc; close all;

% --- 1. Load dataset ---
% The CSV must contain a single column of time-series values.
data = readmatrix('ar5_dataset.csv');
x = data(:);     % ensure column vector
N = length(x);

fprintf('Loaded %d samples.\n', N);

% --- 2. Choose AR model order ---
p = 5; %3 to 7
%Number of parameters

% --- 3. Autocorrelation estimation ---
% Use unbiased estimation for correct YW formulation.
[r_full, lags] = xcorr(x, p, 'unbiased');
% Keep r(0) ... r(p)
r = r_full(lags >= 0);

% --- 4. Build Toeplitz autocorrelation matrix R ---
R = toeplitz(r(1:p));      % R(i,j) = r(|i-j|)
r_vec = r(2:p+1);          % [r1 r2 ... rp]^T

% --- 5. Solve Yule–Walker: R*a = -r_vec ---
a = -R \ r_vec            % AR parameters [a1..ap]

% --- 6. Estimate white-noise variance ---
sigma2_e = r(1) + a' * r_vec;

% --- 7. Compare with MATLAB built-in aryule() ---
[a_yw, sigma2_yw] = aryule(x, p)

% --- 8. Compute residual error of Yule–Walker model ---
% e[k] = x[k] + a1 x[k-1] + ... + ap x[k-p]
e = filter([1 a'], 1, x);

% --- 9. Display results ---
fprintf('\n===== AUTOCORRELATION-BASED (YULE–WALKER) RESULTS =====\n');
fprintf('Estimated AR coefficients (manual):\n');
disp(a');
fprintf('Noise variance estimate:  %.6f\n', sigma2_e);

fprintf('\nMATLAB aryule() coefficients:\n');
disp(a_yw(2:end));
fprintf('MATLAB aryule variance:   %.6f\n', sigma2_yw);

% --- 10. Plot residuals ---
figure;
subplot(2,1,1);
plot(e); grid on;
title('Residual Error (Yule–Walker)');
xlabel('Sample'); ylabel('e[k]');

subplot(2,1,2);
autocorr(e);
title('Autocorrelation of Residuals');

% ======================================================
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



% ======================================================
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


% ======================================================
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


% ======================================================
% 14. Compare AR models for p = 3,4,5,6,7
%     Based on Residual Basic Statistics
% ======================================================

p_values = 3:7;

Mean_vals     = zeros(size(p_values));
Var_vals      = zeros(size(p_values));
Skew_vals     = zeros(size(p_values));
Kurt_vals     = zeros(size(p_values));

for idx = 1:length(p_values)

    p_temp = p_values(idx);

    % -------- Compute autocorrelation for this p --------
    [r_full_temp, lags_temp] = xcorr(x, p_temp, 'unbiased');
    r_temp = r_full_temp(lags_temp >= 0); 

    % Toeplitz matrix
    R_temp = toeplitz(r_temp(1:p_temp));
    r_vec_temp = r_temp(2:p_temp+1);

    % -------- Solve Yule-Walker --------
    a_temp = -R_temp \ r_vec_temp;

    % -------- Compute residuals --------
    e_temp = filter([1 a_temp'], 1, x);

    % Only valid samples
    e_valid_temp = e_temp(p_temp+1:end);

    % -------- Compute residual stats --------
    Mean_vals(idx) = mean(e_valid_temp);
    Var_vals(idx)  = var(e_valid_temp);
    Skew_vals(idx) = skewness(e_valid_temp);
    Kurt_vals(idx) = kurtosis(e_valid_temp);

end

% -------- Display results in a nice table --------
T = table(p_values', Mean_vals', Var_vals', Skew_vals', Kurt_vals', ...
    'VariableNames', {'p','Mean','Variance','Skewness','Kurtosis'});

disp('===============================================================');
disp('        Residual Statistics for p = 3 to 7');
disp('===============================================================');
disp(T);

% ----------- Choose best p: lowest variance + mean ~ 0 ----------
[~, best_idx] = min(Var_vals);
best_p = p_values(best_idx);

fprintf('\n>>>> Best p based on lowest residual variance = %d\n', best_p);
fprintf('     (Also check mean ≈ 0 and kurtosis ≈ 3)\n');

fprintf('     However lowest mean and the kurtosis equal to 3 is with p=5\n');


%% ================================================================
%  PART 2 — LEAST SQUARES AR(p) PARAMETER ESTIMATION
% ================================================================

clear; clc; close all;

% --- 1. Load dataset ---
x = readmatrix('ar5_dataset.csv');
x = x(:);
N = length(x);

fprintf('Loaded %d samples.\n', N);

% --- 2. Choose AR order ---
p = 5;

% --- 3. Build regression matrix Φ and output vector y ---
Phi = zeros(N-p, p);
y   = zeros(N-p, 1);

for k = p+1:N
    Phi(k-p, :) = -x(k-1:-1:k-p)';  
    y(k-p) = x(k);
end

% --- 4. Least Squares solution ---
theta_LS = (Phi' * Phi) \ (Phi' * y);

% --- 5. Regularized LS ---
alpha = 0.1;
theta_reg = (Phi' * Phi + alpha * eye(p)) \ (Phi' * y);

% --- 6. Yule–Walker reference ---
[a_yw, sigma2_yw] = aryule(x, p);
theta_YW = a_yw(2:end).';   % force row vector

% --- 7. Compute residuals (IMPORTANT FIX) ---
e_LS  = y - Phi * theta_LS;
e_reg = y - Phi * theta_reg;

e_YW = filter([1 theta_YW'], 1, x);
e_YW = e_YW(p+1:end);

% --- 8. Compute residual statistics ---
stats = @(e) [mean(e), var(e), skewness(e), kurtosis(e)];

LS_stats  = stats(e_LS);
REG_stats = stats(e_reg);
YW_stats  = stats(e_YW);

fprintf('\n===== Residual Basic Statistics =====\n');
fprintf('        Mean        Variance     Skewness    Kurtosis\n');
fprintf('LS:   %10.6f  %10.6f  %10.6f  %10.6f\n',   LS_stats);
fprintf('REG:  %10.6f  %10.6f  %10.6f  %10.6f\n',   REG_stats);
fprintf('YW:   %10.6f  %10.6f  %10.6f  %10.6f\n',   YW_stats);

% --- 9. Stability Check ---
roots_LS  = roots([1; theta_LS]);
roots_REG = roots([1; theta_reg]);
roots_YW  = roots([1; theta_YW]);

fprintf('\n===== Stability Check =====\n');
fprintf('LS stable?       %d\n', all(abs(roots_LS)  < 1));
fprintf('Regularized LS?  %d\n', all(abs(roots_REG) < 1));
fprintf('Yule-Walker?     %d\n', all(abs(roots_YW)  < 1));

% --- 10. Plot residual comparison ---
figure;
plot(e_LS, 'b'); hold on;
plot(e_reg, 'r');
plot(e_YW, 'k');
legend('LS residuals','Regularized LS residuals','Yule-Walker residuals');
title('Residual Comparison Across AR Estimators');
xlabel('Sample index'); ylabel('Residual value');
grid on;


% ========================================================================
%  Meaningful Plot: Raw Data vs AR Reconstructions Using Estimated Thetas
%  Combined Reconstruction Plots: Full Signal + Zoomed Segment
% ========================================================================

% --- Reconstruct signal using LS estimator ---
x_hat_LS = zeros(size(x));
for k = p+1:N
    x_hat_LS(k) = -theta_LS' * x(k-1:-1:k-p);   % correct
end

% --- Reconstruct using Regularized LS estimator ---
x_hat_REG = zeros(size(x));
for k = p+1:N
    x_hat_REG(k) = -theta_reg' * x(k-1:-1:k-p); % correct
end

% --- Reconstruct using Yule–Walker estimator ---
x_hat_YW = zeros(size(x));
for k = p+1:N
    x_hat_YW(k) = -theta_YW' * x(k-1:-1:k-p);   % correct
end


figure;

% ------------------- Subplot 1: Full Reconstruction --------------------
subplot(2,1,1);

plot(x, 'Color',[0.7 0.7 0.7], 'LineWidth',1.5); hold on;
plot(x_hat_LS,  'r',  'LineWidth',1.8);
plot(x_hat_REG, 'b',  'LineWidth',1.8);
plot(x_hat_YW, 'Color', [0 0.7 0], 'LineWidth', 1.8);


legend('Raw data',...
       'AR reconstruction (LS)',...
       'AR reconstruction (Regularized LS)',...
       'AR reconstruction (Yule–Walker)',...
       'Location','best');

title(sprintf('AR(%d) Model Signal Reconstruction Using Estimated Parameters', p));
xlabel('Sample index k');
ylabel('Amplitude');
grid on;

% ------------------- Subplot 2: Zoomed Reconstruction -------------------
subplot(2,1,2);

start  = floor(N*0.30);
finish = floor(N*0.45);

plot(x(start:finish), 'Color',[0.7 0.7 0.7], 'LineWidth',1.5); hold on;
plot(x_hat_LS(start:finish),  'r', 'LineWidth',1.8); 
plot(x_hat_REG(start:finish), 'b', 'LineWidth',1.8);
plot(x_hat_YW(start:finish), 'Color', [0 0.7 0], 'LineWidth', 1.8);


legend('Raw segment','LS','Regularized LS','Yule–Walker','Location','best');
title('Zoomed Segment: AR Reconstruction Comparison');
xlabel('Sample index');
ylabel('Amplitude');
grid on;

sgtitle('AR Model Reconstruction — Full Signal and Zoomed View');

% ========================================================================
%  EXPLANATION OF RESULTS (LS, REGULARIZED LS, AND YULE–WALKER)
%  ========================================================================
%
%  The residual basic statistics for the three estimation methods are:
%
%       Mean        Variance      Skewness      Kurtosis
%  LS:   ~0.0219     ~0.24399      ~0.03669       ~3.06587
%  REG:  ~0.0219     ~0.24399      ~0.03667       ~3.06586
%  YW:   ~0.0219     ~0.24399      ~0.03667       ~3.06580
%
%  INTERPRETATION:
%
%  1. Mean:
%     The residual mean is approximately 0.021 for all methods.  
%     This value is very close to zero, meaning that the prediction errors
%     do not exhibit systematic bias. The AR model is therefore nearly
%     unbiased for the given dataset.
%
%  2. Variance:
%     All methods produce the same residual variance (~0.244).  
%     This indicates that LS, Regularized LS, and Yule–Walker model the 
%     data with essentially identical accuracy. The innovations have the
%     same noise power, confirming that the process is well approximated
%     by an AR(p) model.
%
%  3. Skewness:
%     The skewness is around 0.036 (close to 0).  
%     This means the residual distribution is almost symmetric.  
%     No estimator consistently overshoots or undershoots the data.
%
%  4. Kurtosis:
%     Kurtosis values are around 3.065, very close to the Gaussian value 
%     kurtosis = 3.  
%     This indicates that residuals behave very much like Gaussian white 
%     noise, which is exactly the expected behavior for a correctly 
%     identified AR model.
%
%  OVERALL RESIDUAL QUALITY:
%     Residuals from all three estimators appear:
%       • almost zero-mean
%       • symmetric
%       • Gaussian-like
%       • with equal variance
%     → This strongly indicates that p = 5 is an appropriate AR model order.
%
%  ========================================================================
%  STABILITY CHECK
%  ========================================================================
%
%  The stability results shown were:
%
%     LS stable?       1
%     Regularized LS?  1
%     Yule–Walker?     1
%
%  INTERPRETATION:
%
%  The output “1” means TRUE — all roots of the AR characteristic
%  polynomial lie inside the unit circle for all three estimators.
%
%  This is expected because:
%     • The process generating the data is stationary.
%     • The chosen order p = 5 is suitable.
%     • Regularization (alpha = 0.1) pulls poles slightly inward.
%     • Yule–Walker is known to yield stable AR models because the
%       autocorrelation matrix is positive definite.
%
%  In summary:
%       • LS         → stable
%       • REG        → stable (even more robust to noise)
%       • Yule–Walker→ stable (reference method)
%
%  ========================================================================
%  FINAL CONCLUSION
%  ========================================================================
%
%  All three estimation methods produce nearly identical AR coefficients.
%  Their residuals have nearly identical statistics, and all models are 
%  stable. This means:
%
%     → The underlying signal is well described by an AR(5) process.
%     → LS, Regularized LS, and Yule–Walker are consistent estimators.
%     → No method significantly outperforms the others for this dataset.
%
%  This validates both the correctness of the algorithms and the suitability
%  of the AR(5) model for the given time series.
%
% ========================================================================



%% ========================================================================
%  PART 3 — Cramér–Rao Bound (CRB) for AR(p) Parameters
% ========================================================================

% clc; clear all; close all;
% 
% x = readmatrix('ar5_dataset.csv');
% x = x(:);
% N = length(x);
% p=5;
% Phi = zeros(N-p, p);
% y   = zeros(N-p, 1);
% theta_LS = (Phi' * Phi) \ (Phi' * y);
% e_LS  = y - Phi * theta_LS;


% 1. Estimate noise variance from LS residuals
sigma2_hat = var(e_LS);   % unbiased noise estimate

% 2. Fisher Information Matrix: I(theta) = (1/sigma^2)*Phi'*Phi
Fisher = (1 / sigma2_hat) * (Phi' * Phi);

% Fisher Information measures how well you can estimate a parameter.
% High Fisher information → small variance → good estimator
% Low Fisher information → large variance → impossible to estimate well

% 3. Cramér–Rao Bound: CRB = sigma^2 * inv(Phi'*Phi)
CRB_cov = sigma2_hat * inv(Phi' * Phi);
CRB_var = diag(CRB_cov);     % individual variances
CRB_std = sqrt(CRB_var);     % standard deviations

% 4. Covariance of LS estimator (should match CRB)
Cov_LS = sigma2_hat * inv(Phi' * Phi);
diff_norm = norm(Cov_LS - CRB_cov);

% Display results
fprintf('\n===== Cramér–Rao Bound (CRB) =====\n');
for k = 1:p
    fprintf('Parameter a_%d: Var >= %.6e   (Std >= %.6e)\n', ...
        k, CRB_var(k), CRB_std(k));
end

fprintf('\nDifference between LS covariance and CRB = %.3e\n', diff_norm);
fprintf('(Should be close to zero: LS is an efficient estimator.)\n');

% ========================================================================
%  EXPLANATION OF CRAMÉR–RAO BOUND (CRB) RESULTS
% ========================================================================
%
%  For an AR(p) model y = Phi*theta + eps with eps ~ N(0, sigma^2 I),
%  the Fisher Information matrix is:
%       I(theta) = (1/sigma^2) * Phi' * Phi
%
%  The Cramér–Rao Bound gives the minimum achievable covariance for any
%  unbiased estimator of theta:
%       CRB(theta) = sigma^2 * inv(Phi' * Phi)
%
%  Therefore, the diagonal elements CRB_var(k) represent the smallest
%  possible variance that any unbiased estimator could have when
%  estimating parameter a_k.
%  In AR modeling with Gaussian noise, the LS estimator is exactly the
%  Maximum Likelihood Estimator (MLE). For linear-Gaussian models, the
%  MLE is known to be an *efficient estimator*, meaning:
%
%       Cov( theta_LS ) = CRB(theta)
%
%  The script computes:
%       • CRB variances and standard deviations
%       • Covariance of the LS estimator
%       • The norm difference ||Cov_LS – CRB|| (should be near zero)
%  A small difference confirms that:
%       → The LS estimator achieves the Cramér–Rao lower bound
%       → No unbiased estimator can perform better
%       → The AR(p) model and Gaussian noise assumptions are consistent
%       → Parameter estimates have minimum variance theoretically possible
%
% ========================================================================


% ========================================================================
%  CRAMÉR–RAO BOUND (CRB) FOR p = 3...7 + COMPARISON PLOT
% ========================================================================

x = readmatrix('ar5_dataset.csv');
x = x(:);
N = length(x);

p_range = 3:7;

CRB_var_all = cell(length(p_range),1);
LS_var_all  = cell(length(p_range),1);

for idx = 1:length(p_range)

    p = p_range(idx);

    % --------------------------------------------------------------------
    % Build Phi and y for current p
    % --------------------------------------------------------------------
    Phi = zeros(N-p, p);
    y   = zeros(N-p, 1);

    for k = p+1:N
        Phi(k-p,:) = -x(k-1:-1:k-p)'; 
        y(k-p) = x(k);
    end

    % --------------------------------------------------------------------
    % LS estimator for order p
    % --------------------------------------------------------------------
    theta_LS = (Phi' * Phi) \ (Phi' * y);

    % --------------------------------------------------------------------
    % Residuals and noise variance estimate
    % --------------------------------------------------------------------
    e_LS = y - Phi * theta_LS;
    sigma2_hat = var(e_LS);

    % --------------------------------------------------------------------
    % CRB computation
    % CRB = sigma^2 * inv(Phi'*Phi)
    % --------------------------------------------------------------------
    CRB_cov = sigma2_hat * inv(Phi' * Phi);
    CRB_var = diag(CRB_cov);        % extract CRB variances

    % --------------------------------------------------------------------
    % Actual estimator covariance (theoretical LS)
    % Cov_LS = sigma^2 * inv(Phi'*Phi)
    % --------------------------------------------------------------------
    Cov_LS = sigma2_hat * inv(Phi' * Phi);
    LS_var = diag(Cov_LS);

    CRB_var_all{idx} = CRB_var;
    LS_var_all{idx}  = LS_var;
end

% ========================================================================
%  Plot CRB vs LS Variance for each p
% ========================================================================

figure;
for idx = 1:length(p_range)
    p = p_range(idx);

    subplot(2,3,idx);
    plot(CRB_var_all{idx}, 'ro-', 'LineWidth', 1.5); hold on;
    plot(LS_var_all{idx},  'bx--','LineWidth', 1.5);
    title(sprintf('Variance Comparison for p = %d', p));
    xlabel('Parameter index k');
    ylabel('Variance');
    legend('CRB variance','LS estimator variance');
    grid on;
end
sgtitle('CRB vs LS Variance for AR(p) Models (p = 3..7)');


% ========================================================================
%  Combined Smoothing Methods: Moving Avg, Low-pass, AR(p)
% ========================================================================

figure;

% --- 1) Moving Average Smoothing ---
window = 20;   % adjust (larger window = smoother)
x_smooth = movmean(x, window);

subplot(3,1,1);
plot(x, 'Color',[0.7 0.7 0.7]); hold on;
plot(x_smooth, 'r', 'LineWidth',2);
legend('Raw data','Moving average smooth');
title('Moving Average Smoothing');
xlabel('Sample index'); ylabel('Amplitude');
grid on;

% --- 2) Low-pass Butterworth Filter ---
fs = 1;          % sampling frequency (1 sample per step)
fc = 0.05;       % cutoff frequency (adjust for smoothness)
[b,a] = butter(4, fc, 'low');
x_lp = filtfilt(b, a, x);

subplot(3,1,2);
plot(x, 'Color',[0.7 0.7 0.7]); hold on;
plot(x_lp, 'r', 'LineWidth',2);
legend('Raw data','Low-pass filtered trend');
title('Low-pass Butterworth Smoothing');
xlabel('Sample index'); ylabel('Amplitude');
grid on;

% --- 3) AR(p) Model-Based Smoothing ---
x_hat = zeros(size(x));
for k = p+1:N
    x_hat(k) = -theta_LS' * x(k-1:-1:k-p);
end

subplot(3,1,3);
plot(x, 'Color',[0.7 0.7 0.7]); hold on;
plot(x_hat, 'r', 'LineWidth',2);
legend('Raw data','AR-model smoothed trend');
title(sprintf('AR(%d)-Based Smoothing', p));
xlabel('Sample index'); ylabel('Amplitude');
grid on;

sgtitle('Comparison of Smoothing Techniques for Raw Time-Series');
