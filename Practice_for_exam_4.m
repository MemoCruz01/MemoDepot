%% ==============================================================
%  EXERCICE 4 – FIRST ORDER PLANT, PI DIGITAL CONTROL
%  G(s) = 2 / (1 + 0.1 s)
%  Desired damping factor λ = 0.42
%  Te = 0.01 s
% ==============================================================
clear; clc; close all;

Te = 0.01;          % Sampling period
s  = tf('s');

%% ==============================================================
% 1) DEFINE THE CONTINUOUS PLANT G(s)
%    G(s) = K / (1 + τ s)  with  K = 2, τ = 0.1
% ==============================================================
K = 2;
tau = 0.1;
G_s = K / (1 + tau*s);

disp('Continuous plant G(s):');
G_s

%% ==============================================================
% 2) DIGITAL TRANSFER FUNCTION OF THE SYSTEM (ZOH)
%    Using the classical ZOH discretization:
%      G(z) = c2d(G(s), Te, 'zoh')
% ==============================================================
G_z = c2d(G_s, Te, 'zoh');
disp('Digital plant G(z):');
G_z

%% ==============================================================
% 3) NICHOLS PLOT OF THE CONTINUOUS PLANT
%    Used to estimate Kp based on damping contour (λ = 0.42)
% ==============================================================
figure;
nichols(G_s); grid on;
title('Nichols Plot of Plant G(s)');

%% ==============================================================
% 4) COMPUTE THE PROPORTIONAL GAIN Kp USING DAMPING SPECIFICATION
%
%  Desired damping factor:
%       λ = 0.42
%
%  Resonance factor:
%       Q = 1 / (2 λ sqrt(1 - λ^2))
%
%  Closed-loop damping specification is satisfied when L(jw)
%  intersects the Q-contour on Nichols.
%
%  We solve analytically:
%
%   Q^2 = |L|^2 / (1 + 2|L| cos φ + |L|^2)
%   with  L(jω) = Kp * G(jω)
%
%   → quadratic equation in Kp:
%
%   (Q^2 − 1)|G|^2 Kp^2 + 2 Q^2 |G| cosφ Kp + Q^2 = 0
%
%  We compute Kp over a frequency grid and pick median.
% ==============================================================
lambda = 0.42;
Q = 1/(2*lambda*sqrt(1-lambda^2));
fprintf("\nDesired damping λ = %.2f → Q = %.4f\n", lambda, Q);

% Frequency grid and frequency response
w = logspace(-1,3,1500);
Gjw = squeeze(freqresp(G_s,w));
magG = abs(Gjw);
phiG = angle(Gjw);

% Solve quadratic at each frequency
Kp_list = NaN(size(w));

for k = 1:length(w)
    g = magG(k);
    c = cos(phiG(k));
    A = (Q^2 - 1) * g^2;
    B = 2 * Q^2 * c * g;
    C = Q^2;
    D = B^2 - 4*A*C;
    if D < 0, continue; end

    K1 = (-B + sqrt(D)) / (2*A);
    K2 = (-B - sqrt(D)) / (2*A);
    K_pos = max([K1 K2]);
    if K_pos > 0
        Kp_list(k) = K_pos;
    end
end

Kp_valid = Kp_list(~isnan(Kp_list));
Kp = median(Kp_valid);

fprintf("Estimated proportional gain Kp = %.4f\n", Kp);

%% ==============================================================
% 5) COMPUTE Ti FROM:
%
%   Ti = 5 / ωr
%
%  For first-order system, the resonant frequency is approximated
%  by the "break frequency" of the system:
%
%     ωr ≈ 1 / τ
%
%  This is valid because first-order systems do NOT have a true
%  resonance peak.
% ==============================================================
omega_r = 1/tau;
Ti = 5 / omega_r;

fprintf("\nResonant frequency ωr = %.4f rad/s\n", omega_r);
fprintf("Integral time Ti = %.4f\n", Ti);

%% ==============================================================
% 6) CONTINUOUS PI CONTROLLER:
%
%       C(s) = Kp ( 1 + 1/(Ti s) )
% ==============================================================
C_s = Kp * (1 + 1/(Ti*s));

disp('Continuous PI controller C(s):');
C_s

%% ==============================================================
% 7) DIGITAL PI CONTROLLER (ZOH)
% ==============================================================
C_z = c2d(C_s, Te, 'zoh');
disp('Digital PI controller C(z):');
C_z

%% --- Robust extraction of digital PI coefficients ---
[numC, denC] = tfdata(C_z,'v');
numC = numC(:).';  % force row

b0 = numC(1);

if length(numC) >= 2
    b1 = numC(2);
else
    b1 = 0;
end

fprintf("\nDigital PI recursive equation:\n");

if b1 == 0
    fprintf("u[k] = u[k-1] + %.6f * e[k]\n", b0);
else
    fprintf("u[k] = u[k-1] + %.6f * e[k] + %.6f * e[k-1]\n", b0, b1);
end

%% ==============================================================
% 8) RECURSIVE EQUATION OF THE CONTROLLER
%
%   If C(z) = (b0 z + b1)/(z - 1),
%   then:
%
%       u[k] = u[k-1] + b0 e[k] + b1 e[k-1]
% ==============================================================
fprintf("\nRecursive PI equation:\n");
fprintf("u[k] = u[k-1] + %.6f e[k] + %.6f e[k-1]\n", b0, b1);

% %% ==============================================================
% % 9) VALIDATION: CLOSED-LOOP RESPONSE (DIGITAL)
% % ==============================================================
% Lz = minreal(C_z * G_z);
% 
% if any(isnan(Lz.Numerator{:})) || any(isnan(Lz.Denominator{:}))
%     error("Open-loop contains NaN. Fix Kp or Ti.");
% end
% 
% CL = feedback(Lz, 1);
% 
% if any(isnan(CL.Numerator{:})) || any(isnan(CL.Denominator{:}))
%     error("Closed-loop contains NaN. Controller not valid.");
% end
% 
% figure; step(CL); grid on;
% figure;
% grid on;
% title('Closed-loop Step Response with Digital PI Controller');