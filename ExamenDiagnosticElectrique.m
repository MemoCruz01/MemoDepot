%% Diagnostic et Observateurs
% -------------------------------------------------------------------------
% System: DC motor (no load torque, no resisting torque)
%
% Electrical equation:
%   u(t) = L di/dt + R i(t) + Ke * w(t)
%
% Mechanical equation:
%   J dw/dt = Kc * i(t)
%
% States chosen in the statement:
%   x1(t) = i(t)   (motor current, A)
%   x2(t) = w(t)   (angular speed, rad/s)
%
% Input:
%   u(t) (supply voltage, V)
%
% Measured outputs:
%   i(t) and w(t), depending on sensor configuration
%
% Lets go for...
%  - Builds state-space model in all sensor configurations
%  - Derives/prints transfer function u -> i
%  - Checks stability numerically
%  - Checks controllability and observability for each C
%  - Simulates response for u = 10V, x(0)=0 on 0..30s
%  - Designs Luenberger observer (current-only measurement), two pole sets
%  - Demonstrates observer performance
%  - Injects sensor biases (faults) and shows estimation error impact
%  - Injects measurement noise and shows estimation error impact
%  - Computes TFs from bias/noise to reconstruction errors (current-only case)

clear; clc; close all;

%% ------------------------------------------------------------------------
% Numerical parameters given in the statement  [oai_citation:2‡Exam_Diag_24-25 (1).pdf](sediment://file_00000000f334720aadb95f9bca45461a)
% -------------------------------------------------------------------------
R  = 0.42;          % Ohm
L  = 0.0027;        % H
Ke = 0.37;          % V/(rad/s)
J  = 2.3e-3;        % kg.m^2
Kc = 1.1;           % N.m/A

fprintf("==== Parameters ====\n");
fprintf("R=%.4f Ohm, L=%.6f H, Ke=%.4f, J=%.6g, Kc=%.4f\n\n", R,L,Ke,J,Kc);

%% ========================================================================
% i. Mise en équation du système : Ecrire le modèle d’état dans chacune des
%    configurations des capteurs suivantes :
%    1) toutes les mesures sont disponibles.
%    2) seule la mesure du courant est disponible.
%    3) seule la mesure de la vitesse angulaire est disponible.
%% ========================================================================

% --- Build A and B from the differential equations ---
% From: u = L di/dt + R i + Ke w
% => di/dt = -(R/L) i - (Ke/L) w + (1/L) u
%
% From: J dw/dt = Kc i
% => dw/dt = (Kc/J) i
%
% Therefore:
A = [ -R/L,    -Ke/L;
       Kc/J,     0   ];

B = [ 1/L;
      0  ];

% Output matrices (C) depend ONLY on which sensors we have.
% Because x = [i; w], measuring i is [1 0]x, measuring w is [0 1]x.
C_all = [1 0;   % y1 = i
         0 1];  % y2 = w

C_i   = [1 0];  % only current measured
C_w   = [0 1];  % only speed measured

D_all = zeros(2,1);
D_i   = 0;
D_w   = 0;

fprintf("==== Q(i) State-space matrices ====\n");
disp("A = "); disp(A);
disp("B = "); disp(B);
disp("C (all measurements) = "); disp(C_all);
disp("C (current only)     = "); disp(C_i);
disp("C (speed only)       = "); disp(C_w);

%% ========================================================================
% ii. Fonction de transfert : En supposant des conditions initiales nulles
%     pour modèle d’état, écrire la fonction de transfert F(s) entre l’entrée
%     u et la sortie i (par calcul littéral puis sur Matlab).
%% ========================================================================

% --- Literal derivation  ---
% Laplace with zero initial conditions:
%   U = (L s + R) I + Ke W
%   J s W = Kc I  =>  W = (Kc/(J s)) I
% Substitute into first:
%   U = (L s + R) I + Ke*(Kc/(J s)) I
% => I/U = 1 / (L s + R + (Ke*Kc)/(J s))
% => I/U = (J s) / (J s (L s + R) + Ke*Kc)
% => F(s) = I(s)/U(s) = (J s) / (LJ s^2 + RJ s + Ke*Kc)

fprintf("\n==== Q(ii) Transfer function F(s)=I(s)/U(s) (literal) ====\n");
fprintf("F(s) = (J*s) / (L*J*s^2 + R*J*s + Ke*Kc)\n");

% --- MATLAB verification using state-space -> transfer function ---
sys_i = ss(A,B,C_i,D_i);
F_tf = tf(sys_i);
fprintf("\n==== Q(ii) Transfer function from MATLAB tf(ss) ====\n");
disp(F_tf);

%% ========================================================================
% iii. Etude de stabilité : Considérant les valeurs numériques, dire si le
%      système est stable ?
%% ========================================================================

fprintf("\n==== Q(iii) Stability check ====\n");
eigA = eig(A);
disp("eig(A) = "); disp(eigA);

if all(real(eigA) < 0)
    fprintf("Conclusion: STABLE (all eigenvalues have negative real part).\n");
else
    fprintf("Conclusion: NOT stable.\n");
end

%% ========================================================================
% iv. Etude de commandabilité et d’observabilité :
%     Dire si ce système est observable dans chacune des trois configurations
%     des capteurs.
%% ========================================================================

fprintf("\n==== Q(iv) Controllability & Observability ====\n");

% Controllability: rank(ctrb) must equal number of states (2)
Rc = rank(ctrb(A,B));
fprintf("rank(ctrb(A,B)) = %d (need 2) -> Controllable? %s\n", Rc, string(Rc==2));

% Observability depends on C
Ro_all = rank(obsv(A,C_all));
Ro_i   = rank(obsv(A,C_i));
Ro_w   = rank(obsv(A,C_w));

fprintf("rank(obsv(A,C_all)) = %d -> Observable (all)? %s\n", Ro_all, string(Ro_all==2));
fprintf("rank(obsv(A,C_i))   = %d -> Observable (i only)? %s\n", Ro_i, string(Ro_i==2));
fprintf("rank(obsv(A,C_w))   = %d -> Observable (w only)? %s\n", Ro_w, string(Ro_w==2));

%% ========================================================================
% v. Simulation du système :
%    Tracer la réponse du système lorsque u = 10V, x(0)=[0;0], horizon 30s.
%% ========================================================================

fprintf("\n==== Q(v) Simulation u=10V, x(0)=[0;0], 0..30s ====\n");
t = 0:0.001:30;                 % fine time grid for smooth plots
u = 10*ones(size(t));           % constant voltage 10V
x0 = [0;0];                     % initial current and speed are zero

sys_all = ss(A,B,C_all,D_all);
[y_true, t_out, x_true] = lsim(sys_all, u, t, x0);

figure;
plot(t_out, y_true(:,1), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('i(t) [A]');
title('Q(v): Current response i(t) for u=10V');

figure;
plot(t_out, y_true(:,2), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('\omega(t) [rad/s]');
title('Q(v): Speed response \omega(t) for u=10V');

%% ========================================================================
% vi. Reconstruction d’état par observateur :
% a) Dans le cas où seul le courant est mesuré, déterminer le gain L
%    de Luenberger pour placement de pôles en {-0.5,-1} puis {-2.5,-3}.
% b) Montrer les performances. xhat(0)=[0.1;0.1]. Différence entre pôles ?
%    Pour la suite: conserver pôles {-2.5,-3}.
%% ========================================================================

fprintf("\n==== Q(vi) Luenberger observer (only current measured) ====\n");

% WHY we can place poles:
% If (A,C_i) is observable -> we can assign eigenvalues of (A - L C_i) arbitrarily.
% MATLAB uses duality:
%   L = place(A', C_i', desired_poles)'

poles_slow = [-0.5 -1.0];
poles_fast = [-2.5 -3.0];

L_slow = place(A', C_i', poles_slow)'; % observer gain for slow poles
L_fast = place(A', C_i', poles_fast)'; % observer gain for fast poles

fprintf("Observer poles %s -> L_slow =\n", mat2str(poles_slow));
disp(L_slow);
fprintf("eig(A - L_slow*C_i) = "); disp(eig(A - L_slow*C_i));

fprintf("\nObserver poles %s -> L_fast =\n", mat2str(poles_fast));
disp(L_fast);
fprintf("eig(A - L_fast*C_i) = "); disp(eig(A - L_fast*C_i));

% --- Performance demonstration ---
% We simulate the real plant and an observer driven by measured current.
% Here, "measured current" is y_i = i(t) (no fault/noise yet).
y_i_meas = y_true(:,1);   % current measurement from simulation

xhat0 = [0.1; 0.1];       % given in statement
[xhat_slow] = simulate_observer_euler(A,B,C_i,L_slow,u,t_out,y_i_meas,xhat0);
[xhat_fast] = simulate_observer_euler(A,B,C_i,L_fast,u,t_out,y_i_meas,xhat0);

% Estimation errors (since x = [i; w])
e_slow = x_true - xhat_slow;
e_fast = x_true - xhat_fast;

figure;
plot(t_out, e_slow(:,1), 'LineWidth', 1.2); hold on;
plot(t_out, e_fast(:,1), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('Error in current  i - î  [A]');
title('Q(vi-b): Current estimation error (slow vs fast poles)');
legend('poles [-0.5 -1]','poles [-2.5 -3]','Location','best');

figure;
plot(t_out, e_slow(:,2), 'LineWidth', 1.2); hold on;
plot(t_out, e_fast(:,2), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('Error in speed  \omega - \omegâ  [rad/s]');
title('Q(vi-b): Speed estimation error (slow vs fast poles)');
legend('poles [-0.5 -1]','poles [-2.5 -3]','Location','best');

fprintf("\nQ(vi-b) Concrete observation:\n");
fprintf("- Fast poles [-2.5 -3] -> faster convergence (errors decay quicker).\n");
fprintf("- Slow poles [-0.5 -1] -> slower convergence (errors decay more slowly).\n");
fprintf("For the rest of the exam we keep poles [-2.5 -3].\n");

% Keep L_fast for the rest
L_keep = L_fast;

%% ========================================================================
% vii. Influence d’un défaut capteur :
% y(t) = Cx(t) + Hm(t)
% a) Add bias:
%    -0.01 A on measured current between 10 and 11s
%    +0.02 rad/s on measured speed between 20 and 21s
%    Observe effect on the observer in each of 3 sensor configurations
%    (recompute observer each time).
% b) Explicit transfer functions from bias -> reconstruction errors of outputs
%    when only current is measured.
%% ========================================================================

fprintf("\n==== Q(vii-a) Sensor bias faults, three sensor configurations ====\n");

% Build bias signals on the full measurement channels [i; w]
m_i = zeros(size(t_out));
m_w = zeros(size(t_out));
m_i(t_out>=10 & t_out<=11) = -0.01;   % bias on current
m_w(t_out>=20 & t_out<=21) = +0.02;   % bias on speed

% Configuration 1: all measurements available y=[i;w]
% Measurement becomes y_meas = [i; w] + [m_i; m_w]
y_meas_all = [y_true(:,1)+m_i, y_true(:,2)+m_w];
C_meas_all = C_all;

% Observer gain for (A,C_all) with the kept poles
L_all = place(A', C_meas_all', poles_fast)';

% Simulate observer driven by 2D measurement
xhat_all = simulate_observer_euler(A,B,C_meas_all,L_all,u,t_out,y_meas_all,[0.1;0.1]);
e_all = x_true - xhat_all;

figure;
plot(t_out, e_all(:,1), 'LineWidth', 1.2); hold on;
plot(t_out, e_all(:,2), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('State estimation error');
title('Q(vii-a): Estimation errors with bias — config 1 (measure i and \omega)');
legend('i - î','\omega - \omegâ','Location','best');

% Configuration 2: only current measured y=i
y_meas_i = y_true(:,1) + m_i;  % only current bias enters
L_i = place(A', C_i', poles_fast)';

xhat_i = simulate_observer_euler(A,B,C_i,L_i,u,t_out,y_meas_i,[0.1;0.1]);
e_i = x_true - xhat_i;

figure;
plot(t_out, e_i(:,1), 'LineWidth', 1.2); hold on;
plot(t_out, e_i(:,2), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('State estimation error');
title('Q(vii-a): Estimation errors with bias — config 2 (measure i only)');
legend('i - î','\omega - \omegâ','Location','best');

% Configuration 3: only speed measured y=w
y_meas_w = y_true(:,2) + m_w;  % only speed bias enters
L_w = place(A', C_w', poles_fast)';

xhat_w = simulate_observer_euler(A,B,C_w,L_w,u,t_out,y_meas_w,[0.1;0.1]);
e_w = x_true - xhat_w;

figure;
plot(t_out, e_w(:,1), 'LineWidth', 1.2); hold on;
plot(t_out, e_w(:,2), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('State estimation error');
title('Q(vii-a): Estimation errors with bias — config 3 (measure \omega only)');
legend('i - î','\omega - \omegâ','Location','best');

fprintf("\nQ(vii-a) Concrete observation:\n");
fprintf("- A bias injected into the measured channel directly perturbs the observer via the innovation (y - yhat).\n");
fprintf("- Config 2 sees only the current bias (10–11s).\n");
fprintf("- Config 3 sees only the speed bias (20–21s).\n");
fprintf("- Config 1 sees both.\n");

%% Q(vii-b) Transfer functions from bias -> reconstruction errors (current-only case)
% When only current is measured:
%   y_meas = C_i x + m(t)  (H=1)
% Observer:
%   xhat_dot = A xhat + B u + L ( y_meas - C_i xhat )
%
% Estimation error x~ = x - xhat satisfies:
%   x~_dot = (A - L C_i) x~ - L * m(t)
%
% Therefore transfer from m -> x~ is:
%   X~(s) = (sI - (A - L C_i))^{-1} * (-L) * M(s)
%
% Because outputs to reconstruct are i and w, and x=[i;w], the reconstruction
% error of outputs equals x~ directly.

fprintf("\n==== Q(vii-b) TFs from bias m -> reconstruction errors (current-only) ====\n");

Acl = A - L_i*C_i;       % error dynamics matrix for current-only observer
Bm  = -L_i;              % m is scalar, so input matrix is (2x1)
Cerr = eye(2);           % output errors are x~ itself (i error and w error)
Dm  = zeros(2,1);

sys_m_to_err = ss(Acl, Bm, Cerr, Dm);
Gm = tf(sys_m_to_err);

disp("Transfer functions from bias m to [i_error; w_error] = ");
disp(Gm);

%% ========================================================================
% viii. Influence du bruit de mesure :
% y(t) = Cx(t) + E b(t) with
%  - b is zero-mean Gaussian, std = 4e-6
%  - E = [0.5  2; 0  1.5]
% a) Plot estimation errors in each sensor configuration.
% b) Explicit TFs from noise -> reconstruction errors (current-only case).
%% ========================================================================

fprintf("\n==== Q(viii) Measurement noise ====\n");

sigma = 4e-6;
E = [0.5 2;
     0   1.5];

% Generate 2 independent noises b1, b2 with given std
rng(1);
b = sigma*randn(numel(t_out),2);  % b(t) is 2D (two noise sources)
noise_y = (E*b')';                % noise_y is Nx2 added to [i;w]

% --- a) errors in each sensor configuration ---

% Config 1 (measure both): y_meas = [i;w] + noise_y
y_meas_all_n = y_true + noise_y;
xhat_all_n = simulate_observer_euler(A,B,C_all,L_all,u,t_out,y_meas_all_n,[0.1;0.1]);
e_all_n = x_true - xhat_all_n;

figure;
plot(t_out, e_all_n(:,1), 'LineWidth', 1.2); hold on;
plot(t_out, e_all_n(:,2), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('Estimation error');
title('Q(viii-a): Estimation errors with noise — config 1 (measure i and \omega)');
legend('i - î','\omega - \omegâ','Location','best');

% Config 2 (measure current only):
% Only the current measurement is available => take first row of noise_y
y_meas_i_n = y_true(:,1) + noise_y(:,1);
xhat_i_n = simulate_observer_euler(A,B,C_i,L_i,u,t_out,y_meas_i_n,[0.1;0.1]);
e_i_n = x_true - xhat_i_n;

figure;
plot(t_out, e_i_n(:,1), 'LineWidth', 1.2); hold on;
plot(t_out, e_i_n(:,2), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('Estimation error');
title('Q(viii-a): Estimation errors with noise — config 2 (measure i only)');
legend('i - î','\omega - \omegâ','Location','best');

% Config 3 (measure speed only):
y_meas_w_n = y_true(:,2) + noise_y(:,2);
xhat_w_n = simulate_observer_euler(A,B,C_w,L_w,u,t_out,y_meas_w_n,[0.1;0.1]);
e_w_n = x_true - xhat_w_n;

figure;
plot(t_out, e_w_n(:,1), 'LineWidth', 1.2); hold on;
plot(t_out, e_w_n(:,2), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('Estimation error');
title('Q(viii-a): Estimation errors with noise — config 3 (measure \omega only)');
legend('i - î','\omega - \omegâ','Location','best');

fprintf("Q(viii-a) Concrete observation:\n");
fprintf("- Noise enters through the innovation term and excites the observer error dynamics.\n");
fprintf("- Faster poles generally amplify noise more (trade-off: speed vs noise robustness).\n");

% --- b) TFs from noise -> reconstruction errors (current-only case) ---
% Current-only measurement: y_meas = C_i x + (E_i * b), where E_i is 1x2
E_i = E(1,:);           % first measurement row
% Error dynamics:
%   x~_dot = (A - L C_i) x~ - L*(E_i b)
% b is 2 inputs -> Bn is 2x2
Bn = -L_i * E_i;        % (2x1)*(1x2) = (2x2)
sys_b_to_err = ss(Acl, Bn, Cerr, zeros(2,2));
Gb = tf(sys_b_to_err);

fprintf("\n==== Q(viii-b) TFs from noise b=[b1;b2] -> reconstruction errors (current-only) ====\n");
disp(Gb);

%% ===================== Local function ===========================
function xhat = simulate_observer_euler(A,B,C,L,u,t,y_meas,xhat0)
% simulate_observer_euler
% WHY Euler?
%  - Extremely transparent: you see directly xhat(k+1)=xhat(k)+dt*xhat_dot
%  - Good for exams and debugging
%
% Observer:
%   xhat_dot = A xhat + B u + L ( y_meas - C xhat )
%
% Works for:
%  - scalar y_meas  (Nx1) when C is 1x2 and L is 2x1
%  - vector y_meas  (Nx2) when C is 2x2 and L is 2x2

dt = t(2)-t(1);
n  = numel(t);

xhat = zeros(n,2);
xhat(1,:) = xhat0(:).';

for k = 1:n-1
    yhat = (C * xhat(k,:).').';          % predicted measurement
    innov = y_meas(k,:) - yhat;          % innovation = y - yhat
    xhat_dot = (A * xhat(k,:).') + (B * u(k)) + (L * innov.');  % observer dynamics
    xhat(k+1,:) = xhat(k,:) + (dt * xhat_dot).';
end
end