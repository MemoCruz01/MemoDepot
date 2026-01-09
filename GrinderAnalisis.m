%% Case study — Grinding-classification process (MATLAB)
% DIAGNOSIS AND OBSERVERS — Case study (Grinding-classification process)
% Numerical values:
%   T1 = 5 min, T2 = 1 min, k1 = 0.5, k2 = 0.3, k3 = 0.1
% Measurement configurations:
%   Case 1: measured overall output flow rate y(t),
%   Case 2: measured output flow rate x1(t) of grinder 1,
%   Case 3: measured output flow rate x2(t) of grinder 2,
%   Case 4: measured overall output flow together with the output flow of grinder 1.
%
% "Whatever the sensor configuration used, the three first outputs will always be reconstructed:
%  the overall flow rate and partial flow rates."  (i.e., [y_overall; x1; x2])
%
% -------------------------------------------------------------------------
clear; clc; close all;

%% Parameters (from statement)
T1 = 5;  T2 = 1;
k1 = 0.5; k2 = 0.3; k3 = 0.1;

% Time base for simulations
t = 0:0.2:30;
dt = t(2)-t(1);

% Input
% The statement says: "Assuming a step input u0 = 1 T/h"
% If flows are in T/min, 1 T/h = 1/60 T/min.
u0_Th = 1;          % 1 T/h (as written)
u0    = u0_Th/60;   % convert to T/min (recommended for consistency)
u     = u0*ones(size(t));

% Initial condition (choose something nonzero to see convergence)
x0 = [1; 1];

%% ========================================================================
% Question 1: System equation
% ▪ Formulate the system operation in symbolic equations.
% ▪ Using as state variables the outlet flow rates of the two grinders,
%   and as the measured variable the system output, write explicitly the
%   state equations of the system in the form:
%   dx(t)/dt = A x(t) + B u(t)
%   y(t)     = C x(t)
%% ========================================================================

% States: x = [x1; x2] (outlet flow of grinder 1 and grinder 2)
% Standard derived model (consistent with the case-study structure):
A = [ -1/T1,   k2/T1;
      k3/T2,  -1/T2 ];

B = [ (1-k1)/T1;
       k1/T2 ];

% Always reconstructed outputs (3 outputs):
% y1 = overall output flow rate
% y2 = x1
% y3 = x2
C_full = [ (1-k3), (1-k2);
           1,      0;
           0,      1 ];
D_full = zeros(3,1);

% Measurement matrices for instrumentation cases:
C1 = [ (1-k3), (1-k2) ];        % Case 1: measure overall output only
C2 = [ 1, 0 ];                  % Case 2: measure x1 only
C3 = [ 0, 1 ];                  % Case 3: measure x2 only
C4 = [ (1-k3), (1-k2);  1, 0 ]; % Case 4: measure [overall; x1]

fprintf("\n=== Q1 ANSWERS ===\n");
disp("A = "); disp(A);
disp("B = "); disp(B);
disp("C_full = [y_overall; x1; x2] = "); disp(C_full);

%% ========================================================================
% Question 2: Stability study
% • Is the system stable?
% • Does stability depend on the numerical values of the parameters ki, and Ti?
% • Check these results on Matlab.
% • In the single-output case (cases 1-3), derive the transfer function:
% • Assuming a step input u0 = 1 T/h, determine the steady-state values.
% • Reproduce and simulate all these results in Matlab; simulate the system’s step.
%% ========================================================================

fprintf("\n=== Q2 ANSWERS: Stability ===\n");
eigA = eig(A);
disp("Eigenvalues(A) = "); disp(eigA);

if all(real(eigA) < 0)
    fprintf("System is STABLE (all eigenvalues have negative real part).\n");
else
    fprintf("System is NOT stable.\n");
end

% Transfer functions for cases 1–3 (SISO)
sys1 = ss(A,B,C1,0);
sys2 = ss(A,B,C2,0);
sys3 = ss(A,B,C3,0);

G1 = tf(sys1); G2 = tf(sys2); G3 = tf(sys3);

fprintf("\n=== Q2 ANSWERS: Transfer functions (single-output cases 1–3) ===\n");
disp("Case 1: G1(s) = Y_overall(s)/U(s) = "); disp(G1);
disp("Case 2: G2(s) = X1(s)/U(s)        = "); disp(G2);
disp("Case 3: G3(s) = X2(s)/U(s)        = "); disp(G3);

% Steady state for step u0 (converted to T/min)
x_ss = -A\(B*u0);
y_ss = C_full*x_ss;

fprintf("\n=== Q2 ANSWERS: Steady-state for step u0 = 1 T/h (converted to 1/60 T/min) ===\n");
disp("x_ss = [x1_ss; x2_ss] = "); disp(x_ss);
disp("[y_overall_ss; x1_ss; x2_ss] = "); disp(y_ss);

% Step simulation (true plant outputs)
sys_full = ss(A,B,C_full,D_full);
[y_true, t_out, x_true] = lsim(sys_full, u, t, x0);

figure;
plot(t_out, y_true, 'LineWidth', 1.2);
grid on; xlabel('t [s]'); ylabel('Output');
title('Q2: Step response (true outputs)'); legend('y_{overall}','x1','x2','Location','best');

%% ========================================================================
% Question 3: Study of controllability and observability
% ▪ Controllability analysis
% • Is the system controllable?
% • Does controllability depend on the numerical values of the parameters ki and Ti?
% ▪ Observability analysis
% • Is the system observable when using only any one of the first three sensors?
% • Does observability depend on the numerical values of the parameters ki and Ti?
%% ========================================================================

fprintf("\n=== Q3 ANSWERS: Controllability ===\n");
Rc = rank(ctrb(A,B));
fprintf("rank(ctrb(A,B)) = %d (need 2)\n", Rc);
fprintf("Controllable? %s\n", string(Rc==2));

fprintf("\n=== Q3 ANSWERS: Observability ===\n");
Ro1 = rank(obsv(A,C1));
Ro2 = rank(obsv(A,C2));
Ro3 = rank(obsv(A,C3));
Ro4 = rank(obsv(A,C4));

fprintf("Case 1 (C1): rank = %d -> Observable? %s\n", Ro1, string(Ro1==2));
fprintf("Case 2 (C2): rank = %d -> Observable? %s\n", Ro2, string(Ro2==2));
fprintf("Case 3 (C3): rank = %d -> Observable? %s\n", Ro3, string(Ro3==2));
fprintf("Case 4 (C4): rank = %d -> Observable? %s\n", Ro4, string(Ro4==2));

%% ========================================================================
% Question 4: State reconstruction via a full-order observer
% ▪ To estimate the full state of the system, use the Luenberger observer with the following structure:
%   d(xhat)/dt = A xhat(t) + B u(t) + L( y(t) - yhat(t) )
% ▪ yhat(t) = C xhat(t)
% ▪ Let x~(t) be the state estimation error defined by: x~(t) = x(t) - xhat(t)
% ▪ Derive the differential equation of the estimation error dynamics.
% ▪ Find how to adjust the gain L so that the estimation error converges to zero.
% ▪ Choose L to assign the observer poles (eigenvalues of A - LC) a priori:
%   {-0.5, -1.0} or {-1.5, -2.0} or {-2.5, -3.0}.
% ▪ In practice, compute L by duality with pole placement on (A', C'):
%   find K s.t. A' - C'K has desired eigenvalues, then L = K'.
%% ========================================================================

% We must pick which measurement is used by the observer.
% Common choice: Case 1 measurement (overall output only)
C_meas = C1;

pole_sets = {[-0.5 -1.0], [-1.5 -2.0], [-2.5 -3.0]};
L_sets = cell(size(pole_sets));

fprintf("\n=== Q4 ANSWERS: Observer gains L for each pole set (using C1) ===\n");
for i=1:numel(pole_sets)
    poles = pole_sets{i};
    L_sets{i} = place(A', C_meas', poles)'; % duality
    fprintf("Poles = %s\n", mat2str(poles));
    disp("L = "); disp(L_sets{i});
    disp("eig(A - L*C) = "); disp(eig(A - L_sets{i}*C_meas));
end

% Choose two pole configurations needed later (Q5/Q6):
L_15_2  = L_sets{2}; % poles [-1.5 -2.0]
L_25_3  = L_sets{3}; % poles [-2.5 -3.0]

%% ========================================================================
% Question 5: Influence of a sensor fault
% We assume sensors 1 and 2 provide biased measurements:
% • +0.5 kg/min between 5.8s and 6.4s for sensor 1,
% • −0.25 kg/min between 11.8s and 12.4s for sensor 2.
% The output equation: y(t) = C x(t) + H m(t), with H = I here.
% Poles placed at {-1.5, -2.0}, then at {-2.5, -3}.
% • Write state estimation error as a function of m(t) and L.
% • Write transfer functions relating m(t) to residual y~(t).
% • Plot reconstruction errors (residuals) for the three outputs.
% Express/plot estimation error and explain how bias can be detected.
%% ========================================================================

% Practical simulation: create faulty measurements for the 3 reconstructed outputs:
%   y_fault = y_true + m(t)  (since H=I)
m = zeros(numel(t),3);
m(t>=5.8 & t<=6.4, 1)  = +0.5;   % sensor 1 bias on y_overall
m(t>=11.8 & t<=12.4,2) = -0.25;  % sensor 2 bias on x1 measurement
% sensor 3: no bias here

y_fault = y_true + m;

% Observer uses ONLY measurement channel Case 1 => y_meas = y_fault(:,1)
y_meas = y_fault(:,1);

% Function: simulate observer for given L (Euler integration for transparency)
simulate_observer = @(Lgain) local_observer_sim(A,B,C_meas,C_full,u,y_meas,t,Lgain);

% Run for the two requested pole sets
[xhat_a, yhat_a, r_a] = simulate_observer(L_15_2);
[xhat_b, yhat_b, r_b] = simulate_observer(L_25_3);

% Residuals: r = y_measured_full - yhat_full
% where y_measured_full is y_fault (all reconstructed outputs are compared)
r_a = y_fault - yhat_a;
r_b = y_fault - yhat_b;

figure;
plot(t, r_a, 'LineWidth', 1.2);
grid on; xlabel('t [s]'); ylabel('Residual');
title('Q5: Residuals with sensor bias (observer poles [-1.5 -2.0])');
legend('r_{overall}','r_{x1}','r_{x2}','Location','best');

figure;
plot(t, r_b, 'LineWidth', 1.2);
grid on; xlabel('t [s]'); ylabel('Residual');
title('Q5: Residuals with sensor bias (observer poles [-2.5 -3.0])');
legend('r_{overall}','r_{x1}','r_{x2}','Location','best');

fprintf("\n=== Q5 CONCRETE ANSWER (how to detect) ===\n");
fprintf("Detect fault by thresholding residuals: |r_i(t)| > gamma for a sustained time.\n");
fprintf("Bias windows should appear clearly in r_overall (5.8–6.4s) and r_x1 (11.8–12.4s).\n");

%% ========================================================================
% Question 6: Influence of measurement noise
% Additive Gaussian noise b(t) with zero mean and variance sigma^2 = 0.0225 on sensors 1,2,3.
% y(t) = Cx(t) + Hm(t) + b(t).
% Poles placed at {-1.5, -2.0}, then at {-2.5, -3}.
% • Write state estimation error as a function of b(t) and L.
% • Write transfer functions relating b(t) to residual y~(t).
% • Plot reconstruction errors (residuals) for the three outputs.
%% ========================================================================

sigma2 = 0.0225;
sigma  = sqrt(sigma2);

rng(1);
b = sigma*randn(numel(t),3);

y_noisy_faulty = y_true + m + b;

y_meas_n = y_noisy_faulty(:,1);

simulate_observer_noise = @(Lgain) local_observer_sim(A,B,C_meas,C_full,u,y_meas_n,t,Lgain);

[xhat_an, yhat_an] = simulate_observer_noise(L_15_2);
[xhat_bn, yhat_bn] = simulate_observer_noise(L_25_3);

r_an = y_noisy_faulty - yhat_an;
r_bn = y_noisy_faulty - yhat_bn;

figure;
plot(t, r_an, 'LineWidth', 1.2);
grid on; xlabel('t [s]'); ylabel('Residual');
title('Q6: Residuals with noise+faults (poles [-1.5 -2.0])');
legend('r_{overall}','r_{x1}','r_{x2}','Location','best');

figure;
plot(t, r_bn, 'LineWidth', 1.2);
grid on; xlabel('t [s]'); ylabel('Residual');
title('Q6: Residuals with noise+faults (poles [-2.5 -3.0])');
legend('r_{overall}','r_{x1}','r_{x2}','Location','best');

fprintf("\n=== Q6 CONCRETE ANSWER ===\n");
fprintf("Noise makes residuals nonzero even without faults -> use thresholds (e.g., 3*sigma) and persistence.\n");
fprintf("Faster poles often increase noise sensitivity but improve fault responsiveness.\n");

%% ========================================================================
% Question 7: Multiple observers and fault signatures
% A bank of observers is used (limit to four observers):
%   - first sensor alone (Case 1),
%   - second sensor alone (Case 2),
%   - first and second sensors (Case 4),
%   - third sensor alone (Case 3).
% Show that this device can detect and locate sensor faults.
%% ========================================================================

fprintf("\n=== Q7: Bank of observers (fault isolation) ===\n");

C_bank = {C1, C2, C4, C3};
names  = {'Obs(C1): overall', 'Obs(C2): x1', 'Obs(C4): [overall;x1]', 'Obs(C3): x2'};

poles_bank = [-1.5 -2.0];
L_bank = cell(size(C_bank));

for i=1:numel(C_bank)
    Ci = C_bank{i};
    try
        L_bank{i} = place(A', Ci', poles_bank)'; % observer gain for that sensor set
        fprintf("%s -> OK, rank(obsv)=%d\n", names{i}, rank(obsv(A,Ci)));
    catch
        fprintf("%s -> NOT designable (not observable with this Ci)\n", names{i});
        L_bank{i} = [];
    end
end

% Simple isolation idea:
% - each observer uses its own measured y_i(t) (with its sensor faults)
% - compare innovation residuals: e_i(t)=y_i(t) - C_i xhat_i(t)
% - the observer that uses the faulty sensor sees a clear innovation signature.

fprintf("Concrete isolation rule:\n");
fprintf("  If sensor 1 is biased: observers using C1 or C4 show large innovation; C2-only and C3-only less affected.\n");
fprintf("  If sensor 2 (x1) is biased: observer C2 and C4 show it strongly.\n");

%% ========================================================================
% Question 8: Influence of a leak
% A leak q(t) may occur at outlet of grinder 2.
% Rewrite as a state equation. Simplify by defining equivalent leak:
%   dq(t)/dt = (-q(t) + f(t))/T2   (as given in statement)
%% ========================================================================

fprintf("\n=== Q8: Augmented state model with leak ===\n");

% Augment with q as an additional state:
% xa = [x1; x2; q], and qdot = (-q + f)/T2, where f(t) unknown (leak forcing).
% One simple physical effect: measured x2_out = x2 - q  (leak subtracts from outlet)
Aa = [A,        [0; 0];   % x dynamics not directly changed here (model choice)
      0 0,     -1/T2];

Ba = [B; 0];              % known input u enters only original x dynamics
Ef = [0; 0; 1/T2];         % unknown input f enters qdot

% Output mapping (example): y_overall uses x2_out = x2 - q
C_full_leak = [ (1-k3), (1-k2), -(1-k2);   % overall depends on x2_out
                1,      0,       0;
                0,      1,      -1];      % x2_out = x2 - q

fprintf("Aa = \n"); disp(Aa);
fprintf("Ba = \n"); disp(Ba);
fprintf("Ef (unknown input) = \n"); disp(Ef);
fprintf("C_full_leak (example) = \n"); disp(C_full_leak);

fprintf("Concrete answer: augment the state with q and add qdot = (-q + f)/T2.\n");

%% ========================================================================
% Question 9: Observer with unknown input
% Leak is an unknown input. Design a state observer with unknown input.
% Show that it is possible to generate an output residual insensitive to the leak
% and sensitive only to sensor/actuator faults.
%% ========================================================================

fprintf("\n=== Q9: Unknown-input observer (UIO) — implementable checklist ===\n");

% A full UIO design requires choosing which outputs are measured (C_meas_leak),
% and checking decoupling rank condition: rank(C_meas_leak*Ef) == rank(Ef).
% If satisfied, one can build a residual that rejects leak input.
%
% Concrete, usable MATLAB checks:
%   rankE = rank(Ef);
%   rankCE = rank(C_meas_leak*Ef);
% If rankCE == rankE -> leak decoupling is feasible in principle.

rankE  = rank(Ef);
fprintf("rank(Ef) = %d\n", rankE);

% Example: suppose we measure [overall; x1] under leak model:
C_meas_leak = C_full_leak(1:2,:); % take first two outputs as measurements (example)
rankCE = rank(C_meas_leak*Ef);
fprintf("Example: rank(C_meas_leak*Ef) = %d\n", rankCE);

if rankCE == rankE
    fprintf("Decoupling condition satisfied (in principle). You can design a UIO to reject the leak.\n");
else
    fprintf("Decoupling condition NOT satisfied with this measurement choice. Change measured outputs.\n");
end

fprintf("Concrete answer: choose measurements so rank(C*Ef)=rank(Ef), then build residual r(t) decoupled from Ef.\n");

%% ---------------- Local function ----------------
function [xhat, yhat_full, innov] = local_observer_sim(A,B,Cmeas,Cfull,u,ymeas,t,L)
    dt = t(2)-t(1);
    n  = length(t);

    xhat = zeros(2,n);
    xhat(:,1) = [0;0];

    innov = zeros(n,1);

    for k=1:n-1
        yhat_meas = Cmeas*xhat(:,k);
        innov(k)  = ymeas(k) - yhat_meas;                 % innovation (y - yhat)
        xhat_dot  = A*xhat(:,k) + B*u(k) + L*innov(k);
        xhat(:,k+1) = xhat(:,k) + dt*xhat_dot;
    end

    yhat_full = (Cfull*xhat)'; % Nx3 reconstructed outputs
end


%% BACKUP NOTES (Concise) — Observability/Controllability + Observer + Stability
% Use this as a quick revision sheet embedded in MATLAB code.
% -------------------------------------------------------------------------
% This file contains ONLY the key exam-style statements and definitions
% (concise), written as comments, so you can copy/paste into any script.

%% 1) Output matrix C and sensor configurations
% - The state-space model is:
%     xdot = A x + B u
%     y    = C x
%
% - IMPORTANT: A and B describe the PHYSICS (plant). They do not change with sensors.
% - Matrix C describes ONLY what is MEASURED (sensor configuration).
%
% - If you define a “bigger C” stacking multiple outputs, e.g.:
%       y_full = [ y_overall; x1; x2 ] = C_full x
%   then different study cases are represented by selecting rows of C_full:
%       Case 1 (measure y_overall only): C = first row of C_full
%       Case 2 (measure x1 only):        C = second row of C_full
%       Case 3 (measure x2 only):        C = third row of C_full
%       Case 4 (measure y_overall & x1): C = [row1; row2]
%
% Concise exam sentence:
%   "Matrix C changes because it represents the measurement configuration.
%    Each study case corresponds to a different set of available sensors,
%    modeled by selecting the appropriate rows of C."

%% 2) Stability (continuous-time LTI)
% - System: xdot = A x + B u
% - Asymptotic stability <=> all eigenvalues of A satisfy Re(lambda) < 0.
%
% Characteristic equation:
%   Eigenvalues lambda are solutions of:
%       det(lambda I - A) = 0
%
% One-line exam summary (for the grinding-classification case):
%   "The eigenvalues of A are always real and strictly negative for Ti > 0
%    and 0 < ki < 1; therefore the system is asymptotically stable,
%    independently of parameter values (within physical ranges)."
%
% MATLAB checks:
%   eig(A)
%   all(real(eig(A)) < 0)

%% 3) Controllability (2-state system)
% Definition:
%   A system is controllable if, using input u(t), we can drive x(t)
%   from any initial state to any final state in finite time.
%
% Test (LTI):
%   Controllability matrix:
%       G = [ B  A*B ]
%   Since n = 2 states:
%       controllable <=> rank(G) = 2
%
% MATLAB checks:
%   G = ctrb(A,B);
%   rank(G)

%% 4) Observability (2-state system)
% Definition:
%   A system is observable if, from measured output y(t),
%   we can reconstruct the internal state x(t).
%
% Test (LTI):
%   Observability matrix:
%       O = [ C
%             C*A ]
%   Since n = 2 states:
%       observable <=> rank(O) = 2
%
% MATLAB checks:
%   O = obsv(A,C);
%   rank(O)
%
% Exam results for the grinding-classification case:
%   - Using only one of the first three sensors (Case 1,2,3): observable (rank=2)
%   - Observability depends on parameters mathematically, but for physical ranges:
%         Ti > 0,  0 < ki < 1
%     the system remains observable.
%
% Compact exam paragraph:
%   "The system is controllable for the given parameters, as the controllability
%    matrix has full rank, although loss of controllability may occur for special
%    parameter combinations. The system is observable using any one of the first
%    three sensors; the corresponding observability matrices have full rank for
%    physically meaningful parameter values."

%% 5) Full-order Luenberger observer
% Observer structure:
%   xhat_dot = A xhat + B u + L ( y - yhat )
%   yhat     = C xhat
%
% Intuition (very short):
%   - A xhat + B u : model prediction (copy of the plant)
%   - (y - yhat)   : innovation (measurement error)
%   - L            : gain correcting the estimate
%   => "A copy of the system, corrected by measurement error."
%
% Estimation error:
%   x_tilde = x - xhat
% Goal: x_tilde(t) -> 0
%
% Error dynamics (key derivation result):
%   x_tilde_dot = (A - L C) x_tilde
%
% Convergence condition:
%   Choose L such that A - L C is stable (its eigenvalues in left half-plane).
%
% Pole placement:
%   If (A,C) is observable, we can assign eigenvalues of (A - L C) arbitrarily.
%
% Practical computation by duality (MATLAB):
%   L = place(A', C', desired_poles)'   % or 'acker' for 2-state
%
% Compact exam paragraph:
%   "A full-order Luenberger observer reconstructs the state.
%    With x~ = x - xhat, the error obeys x~dot = (A - L C)x~.
%    Thus convergence depends only on eigenvalues of (A - L C).
%    Since the system is observable, poles can be assigned by selecting L,
%    computed via dual pole placement on (A^T, C^T)."

%% (Optional) MATLAB one-liners you will repeatedly use
% eig(A)
% all(real(eig(A))<0)
% rank(ctrb(A,B))
% rank(obsv(A,C))
% L = place(A',C',poles)'
% eig(A - L*C)