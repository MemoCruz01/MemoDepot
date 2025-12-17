%%  Continuous to discrete
clc;
lamda = 1.6;
w0= 0.6;
Gs=10  ;

Te=0.01;

%s=tf('s')
num=[3.6 0];
den = [1 1.92 0.36];%[1/(w0)^2 (2*lamda)/w0 1];
Gs = tf(num,den)
Gz = c2d(tf(num,den),Te,'zoh')


clear; clc;

y = [0 0.46 0.30 0.15 0.30 0.61 0.30 0.0 0.0 0.30 0.46 0.15 -0.30 -0.15 0.3 0.15 0.15 0.30 0.61 1.98 6.1 10.06 11.73 12.63 12.84 12.84 12.79 12.76 12.76];
u = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.30 0.30 0.30 0.30 0.30 0.30 0.30 0.30 0.30 0.30 0.30];
t = [0:0.2:5.6];

N=length(y);
Phi= [ -y(2:N-1) -y(1:N-2) -u(2:N-1) -u(1:N-2)];
Y = y(3:N);
theta = [];

theta = Phi\Y%(Phi' * Phi) \ (Phi' * Y);

a1 = theta(1)
a0 = theta(2)
b1 = theta(3)
b0 = theta(4)

%% estimated value 
y_est = [];
y_est(1) = y(1);
y_est(2) = y(2);

for i=3:length(y)
    y_est(i) = -a1* y_est(i-1) - a0* y_est(i-2) + b1* u(i-1) + b0 * u(i-2);
end

figure;
plot(0:length(y)-1, y_est); hold on;
plot(0:length(y)-1, y)

%% Linear Regression Method


LSM = Phi \ y;







