clear; clc; close all;

%% ================================================================
%  PARAMETERS OF THE NONLINEARITY (SATURATION)
% ================================================================
M  = 1;        % saturation output limit
xM = 0.5;      % boundary between linear and saturated zone
k  = M/xM;     % slope inside linear region
% % % 
% % % k=1/0.14
% % % M=k*xM
% % % fprintf("k found qfter looking visuqlly for 0.14 vqlue in x axis = %.3f\n",k);

%% ================================================================
%  LINEAR PLANT  L(s) = K / [ s(1+τ1 s)(1+τ2 s) ]
% ================================================================
K  = 1;
tau1 = 0.1;
tau2 = 1;

s = tf('s');
L = K / ( s*(1+tau1*s)*(1+tau2*s) );

%% ================================================================
%  FREQUENCY RANGE FOR NYQUIST
% ================================================================
w = logspace(-2, 2, 2000);   % 2000 frequency points
Ljw_complex = squeeze(freqresp(L, w));   % complex response
ReL = real(Ljw_complex);
ImL = imag(Ljw_complex);

%% ================================================================
%  AMPLITUDE RANGE FOR DESCRIBING FUNCTION
% ================================================================
A = linspace(0.01, 3, 500);   % amplitudes
N = zeros(size(A));

for i = 1:length(A)
    Ai = A(i);

    if Ai <= xM
        N(i) = k;   % no saturation
    else
        xs = xM/Ai;
        ts = asin(xs);
        N(i) = (2/(pi*Ai)) * ( M*ts + k*xM*sqrt(1-xs^2) );
    end
end

invN = -1 ./ N;   % for Nyquist intersection

%% ================================================================
%  PLOT N(A) and Equivalent Gain 1/|N(A)|
% ================================================================
figure;
subplot(2,1,1);
plot(A, N, 'LineWidth', 2);
grid on;
xlabel('Amplitude A');
ylabel('N(A)');
title('Describing Function N(A) of the Saturation');

subplot(2,1,2);
plot(A, 1./abs(N), 'LineWidth', 2);
grid on;
xlabel('Amplitude A');
ylabel('Equivalent Gain 1/|N(A)|');
title('Equivalent Linear Gain vs Amplitude');

%% ================================================================
%  NYQUIST PLOT + -1/N(A)
% ================================================================
figure; hold on; grid on;

plot(ReL, ImL, 'b', 'LineWidth', 1.5);                % Nyquist curve
plot(invN, zeros(size(invN)), 'r', 'LineWidth', 2);   % -1/N(A)
axis([-3 0.5 -1.5 1.5]);

xlabel('Real axis');
ylabel('Imag axis');
title('Nyquist Plot with -1/N(A) Curve');
legend('Nyquist L(j\omega)', '-1/N(A)');


%No crossing betzeen -1/N(A) zith Nyquist curve assures no osscillation of
%system, therefore stable
%% ================================================================
%  INTERSECTION DETECTION
% ================================================================
ReMat   = ReL(:);          % N×1
ImMat   = ImL(:);          % N×1
invNMat = invN(:).';       % 1×M

distances = abs( (ReMat - invNMat) + 1i * ImMat );

[minD, idx] = min(distances(:));
[row, col] = ind2sub(size(distances), idx);

A_star = A(col);
w_star = w(row);

fprintf("\n===========================================\n");
fprintf(" POSSIBLE LIMIT CYCLE DETECTED\n");
fprintf(" Amplitude  A*  ≈ %.4f\n", A_star);
fprintf(" Frequency ω*  ≈ %.4f rad/s\n", w_star);
fprintf("===========================================\n\n");

% Mark intersection
plot(ReL(row), ImL(row), 'ko', 'MarkerSize', 10, 'LineWidth', 2);
plot(invN(col), 0, 'ks', 'MarkerSize', 10, 'LineWidth', 2);

legend('Nyquist L(j\omega)', '-1/N(A)', ...
       'Intersection on Nyquist', 'Intersection on -1/N(A)');

% %% ================================================================
% %  PLOT LIMIT CYCLE FREQUENCY & AMPLITUDE MARKERS
% % ================================================================
% figure;
% subplot(2,1,1);
% plot(A, N, 'LineWidth', 2); hold on;
% plot(A_star, N(col), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% grid on;
% xlabel('A'); ylabel('N(A)');
% title('Describing Function with Limit Cycle Amplitude Marked');
% 
% subplot(2,1,2);
% plot(w, abs(Ljw_complex), 'b'); hold on;
% plot(w_star, abs(Ljw_complex(row)), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% grid on;
% xlabel('\omega (rad/s)');
% ylabel('|L(j\omega)|');
% title('Linear System Gain Magnitude with Limit Cycle Frequency');

% %% ================================================================
% %  ANIMATION OF -1/N(A) MOVING ALONG REAL AXIS ON NYQUIST PLOT
% % ================================================================
% figure; hold on; grid on;
% 
% % Plot the Nyquist curve once
% plot(ReL, ImL, 'b', 'LineWidth', 1.5);
% axis([-3 0.5 -1.5 1.5]);
% 
% xlabel('Real Axis'); ylabel('Imag Axis');
% title('Animation: Nyquist L(j\omega) and Moving -1/N(A)');
% axis equal;
% 
% % Pre-plot -1/N(A) curve (static)
% plot(invN, zeros(size(invN)), 'r--', 'LineWidth', 1);
% 
% % Prepare animated marker
% h_marker = plot(invN(1), 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% 
% % Optional: create GIF
% makeGIF = true;
% gifname = 'LimitCycle_Animation.gif';
% 
% for i = 1:length(A)
% 
%     % Move marker to -1/N(A(i))
%     set(h_marker, 'XData', invN(i), 'YData', 0);
% 
%     % Highlight intersection if close
%     % Check distance at each frequency
%     dists = abs( (ReL.' - invN(i)) + 1i*ImL.' );
%     [dmin, idx] = min(dists);
% 
%     if dmin < 0.02    % threshold for visual intersection
%         plot(ReL(idx), ImL(idx), 'ko', 'MarkerSize', 10, 'LineWidth', 2);
%     end
% 
%     drawnow;
% 
%     % Save GIF frames
%     if makeGIF
%         frame = getframe(gcf);
%         [imind,cm] = rgb2ind(frame2im(frame),256);
%         if i == 1
%             imwrite(imind,cm,gifname,'gif','DelayTime',0.04,'Loopcount',inf);
%         else
%             imwrite(imind,cm,gifname,'gif','DelayTime',0.04,'WriteMode','append');
%         end
%     end
% end
% 
% disp("Animation completed. GIF saved as " + gifname);

%% ================================================================
%  ADDITIONAL BLOCK: Nyquist using MATLAB's nyquist() function (safe)
% ================================================================
figure; hold on; grid on;

% Nyquist plot using built-in function
[h, wout] = nyquist(L);  % h is complex, wout is frequency vector
h = squeeze(h);           % remove singleton dimensions for SISO

% Plot Nyquist curve
plot(real(h), imag(h), 'b', 'LineWidth', 1.5);

% Overlay -1/N(A) curve
plot(invN, zeros(size(invN)), 'r', 'LineWidth', 2);

% --- Find closest intersection on this Nyquist grid ---
% Compute distance for each frequency in nyquist() output
distances_nyq = abs(real(h) - invN(col) + 1i*imag(h));
[~, idx_nyq] = min(distances_nyq);

% Mark intersection point on nyquist() plot
plot(real(h(idx_nyq)), imag(h(idx_nyq)), 'ko', 'MarkerSize', 10, 'LineWidth', 2);
plot(invN(col), 0, 'ks', 'MarkerSize', 10, 'LineWidth', 2);

xlabel('Real Axis'); ylabel('Imag Axis');
title('Nyquist Plot with -1/N(A) using nyquist() (Safe)');
legend('Nyquist L(j\omega)','-1/N(A)','Intersection on Nyquist','Intersection on -1/N(A)');
axis equal;

