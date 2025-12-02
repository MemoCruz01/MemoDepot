clear; clc;

k = 1:10;
e = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
s = [0, 0.52, 0.9, 1.2, 1.4, 1.55, 1.67, 1.73, 1.82, 1.86];
N = length(s) - 1;
a_values = []; 
b_values = [];

for i=1:9
    y(i) = s(i+1);
    phi(i, 1) = s(i);
    phi(i, 2) = e(i);
end

y = y';

for k=1:N-1
    for j=k+1:N
        PHI = [s(k) e(k); s(j) e(j)];
        Y = [s(k+1) ; s(j+1)];
    
        theta = PHI \ Y;
    
        a_values(end+1) = -theta(1);
        b_values(end+1) = theta(2);
    end

end

fprintf('Values of a = %.4f\n', a_values);
fprintf('Values of b = %.4f\n', b_values);

figure;
plot(1:length(a_values), a_values);
hold on;
plot(1:length(b_values), b_values);

a_est = mean(a_values);
b_est = mean(b_values);

Y_est = [];
Y_est(1) = 0.5;

for k=1:N-1
    Y_est(k+1) = - a_est * Y_est(k) + b_est * e(k);
end

figure;
plot(1:length(Y_est), Y_est); hold on;
plot(1:length(y), y);
grid on;