%% Describing function of a symmetric saturation (first harmonic method)
clear;                                  % Clear all variables from the workspace
close all;                               % Close all open figure windows
clc;                                     % Clear the Command Window

%% Parameters of the sinusoidal input x(t) = A*sin(omega*t)
A = 2;                                   % Input amplitude (this is x1 in your slide)
omega = 5;                               % Input angular frequency [rad/s]
T = 2*pi/omega;                          % Period of the sinusoid

%% Parameters of the saturation nonlinearity
xM = 1;                                  % Saturation threshold on input (breakpoint)
k = 2;                                   % Linear gain (slope) in the non-saturated region
M = k*xM;                                % Output saturation level (w saturates at +/- M)

%% Time vector (use many points per period for good harmonic accuracy)
Npoints = 20000;                         % Number of time samples in one period
t = linspace(0, T, Npoints);             % Time vector covering exactly one period

%% Generate the input signal
x = A*sin(omega*t);                      % Input sinusoid x(t)

%% Apply the saturation nonlinearity to obtain w(t)
w = k*x;                                 % Start with linear behavior w = kx
w(x >  xM) =  M;                         % Saturate to +M when x exceeds +xM
w(x < -xM) = -M;                         % Saturate to -M when x is below -xM

%% Compute the describing function N(A) using the analytical formula
if A <= xM                               % Check if saturation does not occur
    N_formula = k;                       % If A <= xM, the nonlinearity is purely linear
else                                     % If A > xM, saturation occurs
    mu = xM/A;                           % Define mu = xM/A (must be between 0 and 1)
    N_formula = (2*k/pi)*(asin(mu) + mu*sqrt(1 - mu^2)); % Describing function formula
end

%% Compute the describing function numerically from the first harmonic of w(t)
sin_ref = sin(omega*t);                  % Reference sine at the fundamental frequency
a1 = (2/T)*trapz(t, w.*sin_ref);         % First sine Fourier coefficient a1 (numerical integral)
N_numeric = a1/A;                        % Describing function estimate N(A) = a1 / A

%% Display results in the command window
fprintf('Parameters: A = %.4f, xM = %.4f, k = %.4f, M = %.4f\n', A, xM, k, M); % Print parameters
fprintf('N_formula  = %.6f\n', N_formula); % Print analytical N(A)
fprintf('N_numeric  = %.6f\n', N_numeric); % Print numerical N(A)

%% Plot input and output waveforms over one period
figure;                                  % Create a new figure
plot(t, x, 'LineWidth', 1.5);            % Plot input x(t)
hold on;                                 % Keep the plot so we can add more curves
plot(t, w, 'LineWidth', 1.5);            % Plot saturated output w(t)
grid on;                                 % Turn on grid for easier reading
xlabel('Time t [s]');                    % Label x-axis
ylabel('Signal amplitude');              % Label y-axis
title('Input x(t) and saturated output w(t) over one period'); % Plot title
legend('x(t) = A sin(\omega t)', 'w(t) = sat(k x(t))');        % Legend

%% Plot the static nonlinearity curve w = f(x)
figure;                                  % Create a new figure
x_plot = linspace(-1.5*A, 1.5*A, 2000);  % x-axis values for plotting the nonlinearity
w_plot = k*x_plot;                       % Start with linear part
w_plot(x_plot >  xM) =  M;               % Saturate upper part
w_plot(x_plot < -xM) = -M;               % Saturate lower part
plot(x_plot, w_plot, 'LineWidth', 1.5);  % Plot the nonlinearity curve
grid on;                                 % Turn on grid
xlabel('x');                             % Label x-axis
ylabel('w');                             % Label y-axis
title('Saturation nonlinearity w = f(x)'); % Plot title

%% Sweep A to reproduce the describing function curve N(A)
A_values = linspace(0.01, 10*xM, 400);   % Range of amplitudes to evaluate (avoid A=0)
N_values = zeros(size(A_values));        % Preallocate vector for speed

for i = 1:length(A_values)               % Loop over amplitude values
    Ai = A_values(i);                    % Current amplitude
    if Ai <= xM                          % If no saturation occurs
        N_values(i) = k;                 % Equivalent gain is k
    else                                 % If saturation occurs
        mu = xM/Ai;                      % Compute mu = xM/A
        N_values(i) = (2*k/pi)*(asin(mu) + mu*sqrt(1 - mu^2)); % Formula for N(A)
    end
end

figure;                                  % Create a new figure
plot(A_values/xM, N_values/k, 'LineWidth', 1.5); % Normalized plot like in your slide
grid on;                                 % Turn on grid
xlabel('A / x_M');                       % Normalized amplitude axis
ylabel('N(A) / k');                      % Normalized describing function
title('Describing function of symmetric saturation (normalized)'); % Plot title







%% Second Ex
% Describing function of a relay with hysteresis (complex describing function)
clear;                                  % Clear all variables from the workspace
close all;                               % Close all open figure windows
clc;                                     % Clear the Command Window

% Parameters of the sinusoidal input x(t) = A*sin(omega*t)
A = 2;                                   % Input amplitude (x1 in the slide)
omega = 5;                               % Input angular frequency [rad/s]
T = 2*pi/omega;                          % Period of the sinusoid

% Parameters of the relay with hysteresis
M = 1;                                   % Relay output level (w = +/- M)
h = 1;                                   % Total hysteresis width (thresholds are +/- h/2)
th_pos = +h/2;                           % Upper switching threshold
th_neg = -h/2;                           % Lower switching threshold

% Time vector (simulate several periods to reach steady behavior)
Nperiods = 5;                            % Number of periods to simulate
NpointsPerPeriod = 20000;                % Time resolution per period (large for accuracy)
t = linspace(0, Nperiods*T, Nperiods*NpointsPerPeriod); % Time vector

% Generate the input sinusoid
x = A*sin(omega*t);                      % Input signal x(t)

% Simulate the relay with hysteresis (state-machine implementation)
w = zeros(size(t));                      % Preallocate output array w(t)
state = -1;                              % Initial relay state (-1 means w = -M, +1 means w = +M)

for n = 1:length(t)                      % Loop over time samples
    if state == -1                       % If currently at -M
        if x(n) >= th_pos                % Switch to +M when crossing +h/2
            state = +1;                  % Update state
        end
    else                                 % If currently at +M
        if x(n) <= th_neg                % Switch to -M when crossing -h/2
            state = -1;                  % Update state
        end
    end
    w(n) = M*state;                      % Output is +/- M depending on state
end

% Keep only the last period (to avoid transient dependence on initial state)
t_last = t(end-NpointsPerPeriod+1:end);  % Time samples of the final period
x_last = x(end-NpointsPerPeriod+1:end);  % Input during final period
w_last = w(end-NpointsPerPeriod+1:end);  % Output during final period

% Numerically compute the first harmonic of w(t) relative to sin and cos
sin_ref = sin(omega*t_last);             % Sine reference at fundamental frequency
cos_ref = cos(omega*t_last);             % Cosine reference at fundamental frequency

a1 = (2/T)*trapz(t_last, w_last.*cos_ref); % Cosine Fourier coefficient (in-phase with cos)
b1 = (2/T)*trapz(t_last, w_last.*sin_ref); % Sine Fourier coefficient (in-phase with sin)

W1 = a1 - 1j*b1;                         % Complex first-harmonic phasor of w(t)
X1 = -1j*A;                              % Complex phasor of x(t) = A*sin(omega t) (since sin = Im{e^{jωt}})
N_numeric = W1 / X1;                     % Describing function estimate N = W1 / X1

rho_numeric = abs(N_numeric);            % Magnitude of describing function
phi_numeric = angle(N_numeric);          % Phase of describing function (radians)

% Analytical describing function (only valid if A > h/2)
if A <= h/2                              % Check validity condition
    N_formula = NaN;                     % Not defined if no switching happens
    rho_formula = NaN;                   % Not defined
    phi_formula = NaN;                   % Not defined
else
    alpha = asin(h/(2*A));               % Alpha from the slide
    N_formula = (4*M/(pi*A))*exp(-1j*alpha); % Formula from the slide
    rho_formula = abs(N_formula);        % Magnitude from formula
    phi_formula = angle(N_formula);      % Phase from formula
end

% Display numerical vs analytical results
fprintf('A = %.4f, M = %.4f, h = %.4f\n', A, M, h);                % Print parameters
fprintf('N_numeric  = %.6f %+.6fj\n', real(N_numeric), imag(N_numeric)); % Print numerical N
fprintf('rho_numeric = %.6f, phi_numeric = %.6f rad\n', rho_numeric, phi_numeric); % Print magnitude and phase
fprintf('N_formula  = %.6f %+.6fj\n', real(N_formula), imag(N_formula)); % Print formula N
fprintf('rho_formula = %.6f, phi_formula = %.6f rad\n', rho_formula, phi_formula); % Print formula magnitude/phase

% Plot input and relay output over the last period
figure;                                  % Create a new figure
plot(t_last, x_last, 'LineWidth', 1.5);  % Plot input sinusoid
hold on;                                 % Hold plot to add relay output
plot(t_last, w_last, 'LineWidth', 1.5);  % Plot relay output
grid on;                                 % Enable grid
xlabel('Time t [s]');                    % X-axis label
ylabel('Amplitude');                     % Y-axis label
title('Relay with hysteresis: input x(t) and output w(t)'); % Title
legend('x(t)=A sin(\omega t)', 'w(t)=relay+hysteresis');    % Legend

% Sweep A to plot rho(A) and phi(A) like the slide
A_values = linspace(0.01, 10*(h/2), 400); % Range of amplitudes (avoid 0)
rho_values = NaN(size(A_values));         % Preallocate magnitude array
phi_values = NaN(size(A_values));         % Preallocate phase array

for i = 1:length(A_values)                % Loop over amplitudes
    Ai = A_values(i);                     % Current amplitude
    if Ai > h/2                           % Describing function exists only if Ai > h/2
        alpha = asin(h/(2*Ai));           % Alpha from slide
        N_i = (4*M/(pi*Ai))*exp(-1j*alpha); % Describing function at Ai
        rho_values(i) = abs(N_i);         % Magnitude
        phi_values(i) = angle(N_i);       % Phase
    end
end

figure;                                   % Create a new figure
plot(A_values, rho_values, 'LineWidth', 1.5); % Plot magnitude vs A
grid on;                                  % Enable grid
xlabel('A');                              % X-axis label
ylabel('\rho(A) = |N(A)|');               % Y-axis label
title('Magnitude of describing function for relay with hysteresis'); % Title

figure;                                   % Create a new figure
plot(A_values, phi_values, 'LineWidth', 1.5); % Plot phase vs A
grid on;                                  % Enable grid
xlabel('A');                              % X-axis label
ylabel('\phi(A) = \angle N(A) [rad]');    % Y-axis label
title('Phase of describing function for relay with hysteresis'); % Title
