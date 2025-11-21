%% Buro Etude

% ================================================================
%   CONTROL SYSTEM SAMPLING & NUMERICAL CONTROLLER – FULL SOLUTION
%   System: G(p) = 2 / (1 + 0.1 p)
%   Sampling period: Te = 0.01 s
% ================================================================

clear; close all; clc;

% p=tf('s')
% ================================================================
% PART I – ZERO-ORDER HOLD (ZOH) APPROXIMATION
% ================================================================

% 1. Use c2d to compute sampled transfer function G(z)
% --------------------------------------------------------
% G(p) = 2 / (1 + 0.1 p)
% We convert the continuous-time transfer function into a discrete-time
% transfer function using Zero-Order Hold (ZOH).
% The MATLAB function c2d automatically applies:
%      G(z) = Z{ ZOH * G(p) }
% where ZOH = (1 - e^{-sTe}) / s
% --------------------------------------------------------

Te = 0.01;
Gp = tf(2, [0.1 1]);      % continuous-time G(p)
Gz_c2d = c2d(Gp, Te, 'zoh')  

% 2. Recall the transfer function of a Zero-Order Hold (ZOH)
% --------------------------------------------------------
% The ZOH has Laplace transform:
%      B(p) = (1 - e^{-p Te}) / p
% And frequency-domain form:
%      B(jω) = (1 - e^{-jωTe}) / (jω)
% --------------------------------------------------------

% 3. Frequency-domain representation of ZOH: B(jw)
% --------------------------------------------------------
% We compute B(jw) numerically for plotting.
% --------------------------------------------------------

w = logspace(0, 4, 2000);      % frequency range
Bjw = (1 - exp(-1j*w*Te))./(1j*w);

[magG, phaseG] = bode(Gp, w);
magG = squeeze(magG);   
phaseG = squeeze(phaseG);
 
% 4. Compare on Nichols/Bode: G(jw) and G(jw)*B(jw)
% --------------------------------------------------------
% G(jw) = continuous-time frequency response
% G(jw)*B(jw) = frequency response * ZOH
% --------------------------------------------------------

GB = magG .* abs(Bjw);                       % magnitude of G(jw)*B(jw)
phaseGB = phaseG + angle(Bjw)*180/pi;        % phase of G(jw)*B(jw)

figure; 
subplot(2,1,1);
semilogx(w, 20*log10(magG), w, 20*log10(GB));
grid on; title('Bode magnitude: G(j\omega) vs G(j\omega)B(j\omega)');
legend('G(j\omega)', 'G(j\omega) B(j\omega)');

subplot(2,1,2);
semilogx(w, phaseG, w, phaseGB);
grid on; title('Bode phase');
legend('G(j\omega)', 'G(j\omega) B(j\omega)');

% 5. Compute the Z-transform manually (using tables)
% --------------------------------------------------------
% G(p) = 2 / (1 + 0.1p)
% Step-invariant discretization yields:
%      G(z) = (1 - e^{-Te/0.1}) * (2) / (z - e^{-Te/0.1})
% --------------------------------------------------------

alpha = exp(-Te/0.1);
num_manual = 2*(1 - alpha);
den_manual = [1  -alpha];

Gz_manual = tf(num_manual, den_manual, Te)

% 6. Compare manual Z-transform with c2d result
% --------------------------------------------------------
% They should match almost perfectly.
% --------------------------------------------------------

figure;
bode(Gz_c2d, Gz_manual);
legend('c2d result', 'manual Z-transform');
grid on;

%% ================================================================
% PART II – NUMERICAL CONTROLLER
% ================================================================

% 1. Derive recursive equations for integration schemes
% --------------------------------------------------------
% Forward rectangular rule:
%     s(k) = s(k-1) + Te * e(k-1)
%
% Backward rectangular rule:
%     s(k) = s(k-1) + Te * e(k)
%
% Trapezoidal rule:
%     s(k) = s(k-1) + (Te/2)*( e(k) + e(k-1) )
% --------------------------------------------------------

% 2. Corresponding Z-transfer functions
% --------------------------------------------------------
% Using Z-transform:
%
% Forward rule:
%   S(z)/E(z) = Te * z^{-1} / (1 - z^{-1})
%
% Backward rule:
%   S(z)/E(z) = Te / (1 - z^{-1})
%
% Trapezoidal rule:
%   S(z)/E(z) = (Te/2)*(1 + z^{-1}) / (1 - z^{-1})
% --------------------------------------------------------

Sy_forward  = tf(Te, [1 -1], Te) * tf(1, [1 0], Te);     % z^{-1} multiplier
Sy_backward = tf(Te, [1 -1], Te);
Sy_trap     = tf([Te/2  Te/2], [1 -1], Te);

% 3. Numerical differentiator (finite difference)
% --------------------------------------------------------
% Backward difference approximation:
%     de(k)/dt ≈ (e(k) - e(k-1)) / Te
%
% Z-transform:
%     D(z)/E(z) = (1 - z^{-1}) / Te
% --------------------------------------------------------

Dz = tf([1 -1], [Te 0], Te)

% 4. Discrete PI and PID controller transfer functions
% --------------------------------------------------------
% PI(s) = Kp + Ki/s
%
% Using trapezoidal integration for the 1/s:
%    I(z) = (Te/2)*(1 + z^{-1})/(1 - z^{-1})
%
% So:
%    C_PI(z) = Kp + Ki * I(z)
%
% PID:
%    D(s) = Kd * s  → numerical derivative above
% --------------------------------------------------------

Kp = 1; Ki = 2; Kd = 0.1;

I_z = Sy_trap;                   % integral part
D_z = Dz;                        % derivative part

C_PI = Kp + Ki * I_z
C_PID = Kp + Ki * I_z + Kd * D_z

% 5. Recursive equations
% --------------------------------------------------------
% PI:
%   u(k) = u(k-1) + Kp*(e(k)-e(k-1)) + Ki*(Te/2)*( e(k) + e(k-1) )
%
% PID:
%   u(k) = u(k-1)
%        + Kp*(e(k)-e(k-1))
%        + Ki*(Te/2)*(e(k)+e(k-1))
%        + Kd*( e(k) - 2e(k-1) + e(k-2) )/Te
% --------------------------------------------------------

%% ============================================================
% PART III — CONTROLLER DESIGN
% Complete MATLAB solution
% ============================================================

clear; close all; clc;

% === System definition ===
Gp = tf(2, [0.1 1]);     % G(p)
Te = 0.1;                % sampling period

% === ZOH discretization ===
Gz = c2d(Gp, Te, 'zoh');

% === Build G(p)*B(p) (continuous ZOH equivalent) ===
s = tf('s');
B = (1 - exp(-s*Te))/s;
G_ZOH_cont = Gp * B;

% === Nichols plot ===
figure;
nichols(G_ZOH_cont);
title('Nichols plot of G(p)B(p)');
grid on;

% === Damping specification ===
lambda = 0.6;
Q = 1 / (2 * lambda * sqrt(1 - lambda^2));

% === Iso-damping contour ===
phi = linspace(-pi, pi, 2000);
G_iso = Q ./ sqrt(1 - 2*(Q^2).*cos(phi) + Q^2);

figure;
nichols(G_ZOH_cont); hold on;
plot(phi*180/pi, 20*log10(G_iso), 'r', 'LineWidth', 2);
legend('G(p)B(p)', 'Iso-damping contour');
title('Nichols locus with iso-contour');
grid on;

% === Find proportional gain Kp so Nichols curve is tangent ===

% Frequency range
w = logspace(-1, 3, 2000);

K_range = logspace(-3, 3, 300);
bestK = NaN;
mindist = inf;

for K = K_range
    [mag, ph] = bode(K*G_ZOH_cont, w);
    mag = squeeze(mag);
    ph = squeeze(ph);    % in degrees

    % interpolate iso-contour at these phases
    iso_mag = interp1(phi*180/pi, 20*log10(G_iso), ph, 'linear', 'extrap');

    % distance
    d = min(abs(20*log10(mag) - iso_mag));

    if d < mindist
        mindist = d;
        bestK = K;
    end
end

Kp = bestK

% === Closed-loop continuous system with P controller ===
Tclosed = feedback(Kp * Gp, 1);

% === Step response of P-controlled system ===
figure;
step(Tclosed);
title('Closed-loop step response with proportional gain');
grid on;

% === Determine resonance frequency wr ===
[magT, phT, wout] = bode(Tclosed);
magT = squeeze(magT);

[peak, idx] = max(magT);
wr = wout(idx)

% === Natural frequency ===
wn = wr / sqrt(1 - lambda^2)

% === Integral time Ti ===
Ti = 5 / wr

% === Continuous PI controller ===
s = tf('s');
C = Kp * (1 + 1/(Ti*s));

% === Digital PI controller (Tustin) ===
Cd = c2d(C, Te, 'tustin')

% === Closed-loop digital control ===
sys_cl = feedback(Cd * Gz, 1);

figure;
step(sys_cl);
title('Closed-loop digital PI controller');
grid on;

% === Recursive equation ===
% From Cd = (b0 + b1 z^-1) / (1 - a1 z^-1):
[b,a] = tfdata(Cd, 'v');

b0 = b(1);
b1 = b(2);
a1 = a(2);

fprintf("\nDigital PI controller difference equation:\n");
fprintf("u[k] = -a1*u[k-1] + b0*e[k] + b1*e[k-1]\n\n");

