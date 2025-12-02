%% ==============================================================
%  EXERCICE 2 – 2nd ORDER PLANT, PI DESIGN (EXAM VERSION)
%  Includes equations and steps as comments.
% --------------------------------------------------------------
%  Continuous plant (given after identification):
%       G(s) = 10 / (p^2/10^2 + 0.1 p + 1)
%           = 1000 / (s^2 + 10 s + 100)
% ==============================================================

clear; clc; close all;

Te = 40e-3;            % Sampling period for the digital controller
s  = tf('s');

%% ==============================================================
% 1) DEFINE THE CONTINUOUS PLANT G(s)
% --------------------------------------------------------------
% Standard 2nd order form:
%   G(s) = K * wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
% Here:  G(s) = 1000 / (s^2 + 10 s + 100)
%   → wn^2 = 100  → wn = 10 rad/s
%   → 2*zeta*wn = 10 → zeta_plant = 10/(2*10) = 0.5
%   → K = 1000 / wn^2 = 1000 / 100 = 10 (but we’ll treat K inside G)
% ==============================================================

G = 1000/(s^2 + 10*s + 100);

%% ==============================================================
% 2) SPECIFICATIONS: DAMPING FACTOR λ (zeta) = 0.7
% --------------------------------------------------------------
% They give desired closed-loop damping:
%   λ = 0.7
% From that we compute the "resonance factor" Q used in Nichols:
%
%   Q = 1 / (2 λ sqrt(1 - λ^2))
%
% This Q defines an iso-damping contour in the Nichols chart.
% ==============================================================

lambda = 0.7;
Q = 1/(2*lambda*sqrt(1 - lambda^2));
fprintf("Desired damping λ = %.2f → Q ≈ %.4f\n", lambda, Q);

%% ==============================================================
% 3) AUTOMATIC COMPUTATION OF Kp USING THE Q-CONDITION
% --------------------------------------------------------------
% Open-loop:   L(s) = Kp * G(s)
% Frequency response: L(jω) = Kp * G(jω)
% Let:
%   |G| = |G(jω)|      (magnitude)
%   φ   = arg G(jω)    (phase)
%
% Closed-loop resonance factor for a given ω satisfies:
%
%   Q^2 = |L|^2 / (1 + 2|L| cos(φ) + |L|^2),   where |L| = Kp * |G|
%
% This is an equation in Kp. Rearranging gives a quadratic:
%
%   (Q^2 - 1)*|G|^2*Kp^2 + 2 Q^2 |G| cos φ * Kp + Q^2 = 0
%
%   → A Kp^2 + B Kp + C = 0, solve for Kp, choose positive root.
%
% We do this on a grid of ω and then take a representative Kp
% (median of valid solutions).
% ==============================================================

% 3.1 Frequency grid and plant frequency response
w = logspace(-1, 3, 2000);          % ω from 0.1 to 1000 rad/s
[Gjw, ~] = freqresp(G, w);
Gjw = squeeze(Gjw);                 % column vector

magG = abs(Gjw);                    % |G(jω)|
phiG = angle(Gjw);                  % ∠G(jω) [rad]

Kp_candidates = NaN(size(w));       % to store possible Kp(ω)

for k = 1:length(w)
    g = magG(k);                    % |G|
    c = cos(phiG(k));               % cos φ

    % Coefficients of A Kp^2 + B Kp + C = 0
    A = (Q^2 - 1) * g^2;
    B = 2 * Q^2 * c * g;
    C = Q^2;

    % Discriminant
    D = B^2 - 4*A*C;
    if D < 0
        continue;                   % no real Kp here
    end

    % Two roots; we only keep positive meaningful Kp
    K1 = (-B + sqrt(D)) / (2*A);
    K2 = (-B - sqrt(D)) / (2*A);
    K_pos = [K1 K2];
    K_pos = K_pos(K_pos > 0);

    if ~isempty(K_pos)
        Kp_candidates(k) = max(K_pos);  % choose larger positive root
    end
end

% Remove NaNs and choose a representative Kp (median of solutions)
Kp_valid = Kp_candidates(~isnan(Kp_candidates));
Kp = median(Kp_valid);

fprintf("\nComputed proportional gain from Q-condition:\n");
fprintf("Kp ≈ %.4f\n", Kp);

%% ==============================================================
% 4) COMPUTE INTEGRAL TIME Ti = 5 / ωr
% --------------------------------------------------------------
% For a 2nd order closed loop:
%   T(s) = ω_n^2 / (s^2 + 2 λ ω_n s + ω_n^2)
%
% Resonant frequency (if λ < 1/√2):
%   ω_r = ω_n * sqrt(1 - 2 λ^2)
%
% Then they impose:  Ti = 5 / ω_r
% ==============================================================

wn = 10;                                      % from plant denominator
omega_r = wn * sqrt(1 - 2*lambda^2);
Ti = 5 / omega_r;

fprintf("\nResonant frequency ωr ≈ %.4f rad/s\n", omega_r);
fprintf("Integral time      Ti ≈ %.4f s\n", Ti);

%% ==============================================================
% 5) CONTINUOUS PI CONTROLLER C(s)
% --------------------------------------------------------------
% Standard PI form:
%   C(s) = Kp ( 1 + 1/(Ti s) )
% ==============================================================

C_s = Kp * (1 + 1/(Ti*s));
disp('Continuous-time PI controller C(s):');
C_s

%% ==============================================================
% 6) DIGITAL PI CONTROLLER C(z) (ZOH DISCRETIZATION)
% --------------------------------------------------------------
% Use c2d with ZOH (as in lecture):
%   C(z) = (b0 z + b1) / (z - 1)
% which corresponds to the difference equation:
%
%   u[k] = u[k-1] + b0 e[k] + b1 e[k-1]
% ==============================================================

C_z = c2d(C_s, Te, 'zoh');
disp('Digital PI controller C(z):');
C_z

[numC, denC] = tfdata(C_z, 'v');
b0 = numC(1);
b1 = numC(2);

fprintf('\nRecursive controller equation:\n');
fprintf('u[k] = u[k-1] + (%.6f) e[k] + (%.6f) e[k-1]\n', b0, b1);

%% ==============================================================
% 7) GRAPH 1 – NICHOLS PLOT WITH DIFFERENT Kp
% --------------------------------------------------------------
% This shows clearly how increasing Kp shifts the Nichols curve
% vertically (in dB) while the phase stays the same.
% ==============================================================

figure; hold on;
K_try = [0.5 1 2 5 10 20];          % some trial gains

for Ki = K_try
    nichols(Ki*G);
end

%ncgrid;
grid on;
title('Effect of K_p on Nichols plot (vertical shift)');
legend('0.5','1','2','5','10','20','Location','best');

%% ==============================================================
% 8) GRAPH 2 – NICHOLS PLOT FOR FINAL DESIGNED Kp
% ==============================================================

figure;
nichols(Kp*G); %ncgrid;
grid on;
title(sprintf('Nichols of L(s) = Kp G(s),  Kp = %.2f', Kp));

%% ==============================================================
% 9) CLOSED-LOOP VALIDATION WITH DIGITAL PI
% --------------------------------------------------------------
% Discretize plant with ZOH for implementation on the computer:
%   G_z(s) = c2d(G(s), Te, 'zoh')
% Closed-loop:
%   T(z) = C(z)G_z(z) / (1 + C(z)G_z(z))
% ==============================================================

G_z = c2d(G, Te, 'zoh');
CL_z = feedback(C_z * G_z, 1);

figure;
step(CL_z);
grid on;
title('Closed-loop step response with designed digital PI');

% (Optional) continuous closed loop for comparison:
CL_s = feedback(C_s * G, 1);
figure;
step(CL_s);
grid on;
title('Continuous closed-loop step response (reference)');