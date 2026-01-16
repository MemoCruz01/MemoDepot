
%% ClosedLoopSaturation
%% Limit cycle prediction for a saturation nonlinearity using describing function + Nyquist
clear;                                      % Clear workspace variables
close all;                                   % Close all figures
clc;                                         % Clear command window

%% -------------------- USER PARAMETERS (EDIT THESE) --------------------
M  = 1;                                      % Saturation output level (+/- M)
xM = 1;                                      % Input breakpoint (saturation starts when |x| >= xM)
k  = M/xM;                                   % Linear slope in the unsaturated region (so that output hits M at xM)

% Define the linear block L(p) (EDIT THIS to match your course example)
s = tf('s');                                 % Define Laplace variable s
L = 5/((s+1)*(s+2));                          % Example stable plant (choose any stable L(s) you want)

% Simulation settings
tEnd = 50;                                   % Simulation duration in seconds
x0plant = [];                                 % Initial plant state (leave empty to use zeros)

%% -------------------- DESCRIBING FUNCTION OF SATURATION --------------------
% For odd symmetric saturation, describing function N(A) is real (no phase).
% It is:
%   N(A)=k                           for A <= xM  (no saturation)
%   N(A)=(2k/pi)*(asin(mu)+mu*sqrt(1-mu^2))  for A > xM, mu=xM/A

Nsat = @(A) (A<=xM).*k + ...                 % Piece for A <= xM (pure linear gain)
            (A>xM).*((2*k/pi).*(asin(xM./A) + (xM./A).*sqrt(1-(xM./A).^2))); % Piece for A > xM

%% -------------------- BUILD THE CRITICAL LOCUS -1/N(A) --------------------
Agrid = linspace(0.05, 10*xM, 600);          % Amplitude grid (avoid 0 to prevent division issues)
Ngrid = arrayfun(Nsat, Agrid);               % Evaluate N(A) on the grid
crit  = -1./Ngrid;                           % Critical locus points on real axis: -1/N(A)

%% -------------------- COMPUTE NYQUIST OF L(jw) --------------------
w = logspace(-3, 3, 4000);                   % Frequency grid (rad/s), wide range
Ljw = squeeze(freqresp(L, w));               % Complex frequency response L(jw) at each w
ReL = real(Ljw);                              % Real part of L(jw)
ImL = imag(Ljw);                              % Imag part of L(jw)

%% -------------------- FIND A LIMIT CYCLE INTERSECTION --------------------
% Condition for intersection with saturation describing function:
%   L(jw0) = -1/N(A)
% Since -1/N(A) is on the negative real axis, we need:
%   Im{L(jw0)} = 0    AND    Re{L(jw0)} < 0
% Then we match magnitudes: N(A) = -1/Re{L(jw0)}

% Find indices where Im{L(jw)} changes sign (potential real-axis crossings)
signChange = find(ImL(1:end-1).*ImL(2:end) <= 0);  % Candidate crossings where Im crosses zero

wCandidates = [];                             % Store candidate frequencies
ReCandidates = [];                            % Store corresponding real values

for idx = signChange.'                        % Loop over candidate crossing segments
    w1 = w(idx);                              % Lower frequency of segment
    w2 = w(idx+1);                            % Upper frequency of segment
    Im1 = ImL(idx);                           % Imag part at w1
    Im2 = ImL(idx+1);                         % Imag part at w2

    if (Im1 == Im2)                            % Avoid degenerate case
        continue;                              % Skip if no variation
    end

    % Linear interpolation to approximate where Im{L(jw)} = 0
    wc = w1 + (w2-w1)*(0 - Im1)/(Im2 - Im1);   % Approx crossing frequency
    Lc = squeeze(freqresp(L, wc));             % L(jwc) value at that frequency (complex)
    if real(Lc) < 0                            % Keep only negative real-axis crossings
        wCandidates(end+1) = wc;               % Store candidate frequency
        ReCandidates(end+1) = real(Lc);        % Store candidate real value
    end
end

if isempty(wCandidates)                        % Check if we found any negative real-axis crossings
    error('No negative real-axis crossing found in Nyquist of L(jw). Try changing L(s) or the frequency range.');
end

% Choose the crossing that is "most negative" (often the relevant one)
[~, bestIdx] = min(ReCandidates);              % Find index of most negative real part
w0 = wCandidates(bestIdx);                     % Selected oscillation frequency candidate
ReAtw0 = ReCandidates(bestIdx);                % Real value of L(jw0) (negative)

% Now solve N(A) = -1 / Re{L(jw0)}
targetN = -1/ReAtw0;                           % Required describing function value

% Define equation f(A)=Nsat(A)-targetN
fA = @(A) Nsat(A) - targetN;                   % Root when N(A)=targetN

% Pick an initial guess for A (start slightly above xM so saturation is active)
Aguess = max(1.2*xM, 0.1);                     % Initial guess for amplitude

% Use fzero to solve for A (limit cycle amplitude x10)
x10 = fzero(fA, Aguess);                       % Solve for amplitude giving intersection

%% -------------------- REPORT THE PREDICTED LIMIT CYCLE --------------------
fprintf('Predicted limit cycle:\n');           % Print header
fprintf('  omega0 = %.6f rad/s\n', w0);         % Print predicted frequency
fprintf('  x10    = %.6f (amplitude)\n', x10);  % Print predicted amplitude
fprintf('  N(x10) = %.6f\n', Nsat(x10));        % Print describing function at that amplitude
fprintf('  L(jw0) = %.6f %+.6fj\n', real(squeeze(freqresp(L,w0))), imag(squeeze(freqresp(L,w0)))); % Print L(jw0)

%% -------------------- PLOT NYQUIST + CRITICAL LOCUS --------------------
figure;                                        % New figure
plot(ReL, ImL, 'LineWidth', 1.5);              % Plot Nyquist of L(jw) for w>=0
hold on;                                       % Hold plot for overlay
plot(crit, zeros(size(crit)), 'r--', 'LineWidth', 1.5); % Plot critical locus -1/N(A) on real axis
plot(-1/Nsat(x10), 0, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Mark the intersection point
grid on;                                       % Enable grid
xlabel('Re\{L(j\omega)\}');                    % Label x-axis
ylabel('Im\{L(j\omega)\}');                    % Label y-axis
title('Nyquist of L(j\omega) and critical locus -1/N(A)'); % Title
legend('Nyquist L(j\omega)', '-1/N(A)', 'Intersection', 'Location', 'Best'); % Legend

%% -------------------- TIME-DOMAIN SIMULATION OF THE NONLINEAR CLOSED LOOP --------------------
% We simulate: w = sat(k*e), e = -y (assuming reference=0), plant: y = L(s)*w
% Convert L(s) to state-space: xdot = Ax + Bw, y = Cx + Dw

sys = ss(L);                                   % Convert transfer function L(s) to state space
Aplant = sys.A;                                % Plant A matrix
Bplant = sys.B;                                % Plant B matrix
Cplant = sys.C;                                % Plant C matrix
Dplant = sys.D;                                % Plant D matrix

if isempty(x0plant)                            % If user didn't specify initial state
    x0plant = zeros(size(Aplant,1),1);         % Use zero initial state
end

% Define the saturation function: sat(v)=clip(k*v, +/-M)
satfun = @(v) max(min(k*v, M), -M);             % Saturation with slope k and limits +/-M

% ODE for closed-loop dynamics
odefun = @(t,x) (Aplant*x + Bplant*satfun(- (Cplant*x + Dplant*0))); % Placeholder, we will compute w properly below

% We implement properly by computing y=Cx + D*w, but w depends on y -> algebraic loop if D!=0.
% To avoid issues, we assume D=0 (common in proper plants). If D != 0, you need a small delay or solve algebraically.
if any(abs(Dplant) > 1e-12)                    % Check if D matrix is nonzero
    warning('D is nonzero; simulation may require algebraic loop handling. Consider using a proper plant with D=0.'); % Warn user
end

% Define the ODE with w computed from y (assuming D=0 for simplicity)
odefun = @(t,x) (Aplant*x + Bplant*satfun(-(Cplant*x))); % xdot = Ax + B*w, w = sat(-y)

% Integrate the ODE
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);   % ODE solver tolerances for accuracy
[tSim, xSim] = ode45(odefun, [0 tEnd], x0plant, opts); % Solve ODE over [0, tEnd]

% Compute output y(t) and input w(t)
ySim = (Cplant*xSim.').';                      % Compute y(t)=C x(t)
wSim = satfun(-ySim);                          % Compute w(t)=sat(-y(t))

%% -------------------- PLOT TIME RESPONSE AND ESTIMATE OSCILLATION --------------------
figure;                                        % New figure
plot(tSim, ySim, 'LineWidth', 1.2);            % Plot y(t)
grid on;                                       % Grid on
xlabel('Time t [s]');                          % X label
ylabel('y(t)');                                 % Y label
title('Nonlinear closed-loop response (should converge to a limit cycle)'); % Title

% Estimate steady-state amplitude from the last portion of the signal
idxSS = tSim > (0.7*tEnd);                      % Indices for last 30% of simulation
ySS = ySim(idxSS);                              % Steady-state segment
ampEst = 0.5*(max(ySS)-min(ySS));               % Amplitude estimate = (peak-to-peak)/2

fprintf('Time-domain estimated amplitude ~ %.6f\n', ampEst); % Print estimated amplitude
fprintf('Describing-function predicted amplitude x10 = %.6f\n', x10); % Print predicted amplitude
