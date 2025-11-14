%% Excercise 4
close all; clear all; clc;



%sys=tf(num,den)

p=tf('s');
G=2/(0.1*p+1);  
% figure;
% bode(G)


%% Shanon
wmax=100
fmax = wmax / (2*pi);
%Theo Shannon
fe_shannon=2*fmax
Te_shannon=1/fe_shannon % Periode a 31ms

%On pose  Te=40ms
fe=1/(40e-3)
fmax=fe/2
wmax=fmax*2*pi

%%
Te=40e-3;
Gc=c2d(G,Te);
figure;
bode(G);
hold all;
bode(Gc);
Te=3e-3;
Gc=c2d(G,Te);
hold all
bode(Gc);
legend('Cont','Discrete 40ms','Discrete 3ms')

%%

figure; 
nichols(G)

hold all
nichols(Gc)
Kp=10^((30-4.5)/20);
nichols(Kp*Gc)
ngrid
legend('Cont','Discrete','Kp*Gz')


