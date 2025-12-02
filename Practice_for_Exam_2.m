clear; clc;

s = [0 0.2 0.55 0.77 0.74];
e = [1 0.8 0.45 0.23 0.26];

f = [s(1) 0 e(1) 0;
    s(2) s(1) e(2) e(1);
    s(3) s(2) e(3) e(2);
    s(4) s(3) e(4) e(3);
    ]

y = [s(2); s(3); s(4); s(5)]

theta = [];

theta = inv(f) * y

a1 = -theta(1);
a0 = -theta(2);
b1 = theta(3);
b0 = theta(4);

%%estimated value 
y_est = [];
y_est(1) = 0.01;
y_est(2) = s(2);

for i=3:length(s)
    y_est(i) = -a1* y_est(i-1) - a0* y_est(i-2) + b1* e(i-1) + b0 * e(i-2);
end

figure;
plot(0:length(s)-1, y_est); hold on;
plot(0:length(s)-1, s)

%% Linear Regression Method

phi = f;

LSM = phi \ y;



