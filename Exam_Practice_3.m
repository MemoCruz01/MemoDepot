clear; clc; close all;

%% ============================================================
% QUESTION 1 and 2
% ============================================================


% Ts = 0.09;
tau = 0.1;

num=2;
den=[0.1 1];
G = tf(num,den) 
% c2d(G,Ts,'zoh')
% 
% figure;
% nichols(G); grid on;

%% ============================================================
% PRACTICAL METHOD: Try several Kp values
% and visually find which intersects desired damping contour.
% ============================================================

Kp_values = [0.1 0.2 0.5 1 2 5 10];

figure; hold on;
for Kp = Kp_values
    L = Kp * G;          % OPEN-LOOP L(s) = Kp * G(s)
    nichols(L);          % plot for each Kp
end

% ncgrid;                  % overlays closed-loop damping/overshoot contours
grid on;
title('Nichols plot of L(s)=Kp·G(s) for different Kp values');
legend('Kp=0.1','Kp=0.2','Kp=0.5','Kp=1','Kp=2','Kp=5','Kp=10', ...
       'Location','best');


%% ============================================================
% QUESTION 3 — Estimate Kp for desired damping λ = 0.42
% Using corner frequency instead of ωr (first-order system)
% ============================================================

% --- Plant definition ---
K   = 1;
s = tf('s');
G = K / (1 + tau*s);

% --- Desired damping ---
lambda = 0.42;
Q = 1 / (2*lambda*sqrt(1 - lambda^2));   % resonance factor (conceptual)
fprintf('Desired damping λ = %.2f  →  Q ≈ %.3f\n', lambda, Q);

% NOTE:
% For a first-order system, there is NO resonant frequency.
% We therefore use the CORNER FREQUENCY:
wc = 1/tau;
fprintf('Corner frequency (1/tau): ωc = %.2f rad/s\n', wc);

% ============================================================
% Evaluate plant at the CORNER frequency
% ============================================================
[mag_c, phase_c] = bode(G, wc);
mag_c   = squeeze(mag_c);      % magnitude at corner frequency
phase_c = squeeze(phase_c);    % phase in degrees

fprintf('At ωc = %.2f rad/s:  |G(jωc)| = %.4f   phase = %.2f°\n', ...
        wc, mag_c, phase_c);


%% ============================================================
% Compute the OPEN-LOOP magnitude needed at ωc to satisfy Q
% (THEORETICAL APPROXIMATION)
% ============================================================

% Q-condition applied at CORNER frequency:
% Q^2 = |L|^2 / (1 + 2|L|cos(phase) + |L|^2)
% Solve for required |L| (target magnitude)

phi_rad = phase_c * pi/180;      % convert to radians
Q2 = Q^2;

% Solve quadratic equation for |L|
% Let Lmag = |L(jωc)| needed
% (Q2 - 1)*Lmag^2 + 2*Q2*cos(phi)*Lmag + Q2 = 0

a = (Q2 - 1);
b = 2*Q2*cos(phi_rad);
c = Q2;

disc = b^2 - 4*a*c;
%Solving second order equation
if disc < 0
    disp('No real solution for required |L| (expected for first-order systems).');
    L_required = NaN;
else
    Lsol = (-b + sqrt(disc)) / (2*a);
    Lalt = (-b - sqrt(disc)) / (2*a);
    L_required = max([Lsol Lalt]);   % positive solution
end

fprintf('Required |L(jωc)| ≈ %.4f  (if valid)\n', L_required);


%% ============================================================
% Estimate Kp from required magnitude (if meaningful)
% ============================================================

if isnan(L_required)
    fprintf('No valid analytic Kp — choose Kp visually from Nichols plot.\n');
else
    Kp_est = L_required / mag_c;
    fprintf('Estimated proportional gain Kp (analytic) = %.4f\n', Kp_est);
end




%% ============================================================
%  EXERCISE 4 — Questions 5, 6, 7 and 8
% ============================================================

Ti=5/wc;

%% =============================
% Known values (INSERT YOURS)
% =============================
Te = 0.01;      % Sampling period
Kp = 1;         % <-- Replace by your value from Question 3

s = tf('s');

%% ============================================================
%  QUESTION 5 – Continuous PI Controller C(p)
% ============================================================

% Continuous-time PI:
% C(p) = Kp ( 1 + 1/(Ti * p) )
C_p = Kp * (1 + 1/(Ti*s));

disp('Continuous PI controller Cp(s):')
C_p


%% ============================================================
%  QUESTION 6 – Digital PI Controller C(z) (ZOH)
% ============================================================

% Discretize the continuous PI controller using ZOH:
C_z = c2d(C_p, Te, 'zoh');

disp('Digital PI controller C(z):')
C_z


%% ============================================================
%  QUESTION 7 – Recursive (Difference) Equation
% ============================================================

% Extract numerator and denominator
[numC, denC] = tfdata(C_z, 'v');

% For a PI discretized with ZOH:
% C(z) = (b0 z + b1) / (z - 1)
% → u[k] = u[k-1] + b0 e[k] + b1 e[k-1]

b0 = numC(1);
b1 = numC(2);

fprintf('\nRecursive PI controller equation:\n');
fprintf('u[k] = u[k-1] + (%.6f) e[k] + (%.6f) e[k-1]\n', b0, b1);


%% ============================================================
%  QUESTION 8 – Closed-loop Validation (Discrete Plant)
% ============================================================

% Continuous plant:
% G(s) = 2 / (1 + 0.1*s)
G_s = 2 / (1 + 0.1*s);

% Discretize plant with ZOH:
G_z = c2d(G_s, Te, 'zoh');

% Closed-loop transfer function:
% T(z) = (C(z) G(z)) / (1 + C(z) G(z))
CL = feedback(C_z * G_z, 1);

% Step response plot:
figure;
step(CL);
grid on;
title('Closed-loop Step Response with Digital PI Controller');


%% ============================================================
%  OPTIONAL: Compare open-loop response
% ============================================================

figure;
bode(C_z * G_z); grid on;
title('Bode Plot of Open-loop (C(z) * G(z))');


%% ============================================================
%  OPTIONAL: Manual simulation of PI recursion + plant
% ============================================================

% %% First-order plant recursion
% [numG, denG] = tfdata(G_z, 'v');
% b0p = numG(1);
% b1p = numG(2);
% a1p = denG(2);
% 
% for k = 2:6
%     % Error
%     e(k) = r(k) - y(k);
% 
%     % Controller recursion
%     u(k) = u(k-1) + b0*e(k) + b1*e(k-1);
% 
%     % Plant recursion (first-order)
%     y(k) = -a1p * y(k-1) + b0p*u(k) + b1p*u(k-1);
% end
% figure;
% plot(y, 'LineWidth', 1.8); grid on;
% title('Manual Simulation of Closed-loop Response');
% xlabel('k'); ylabel('y[k]');