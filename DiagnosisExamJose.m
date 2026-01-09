% Master MARS 2025-2026
% Jose
% Exam 09/01/26

clc; clear all; close all;
% System parameters
k1=0.1;    % in N/m
k2=0.3;    % in N/m
alpha=1.5; % in Ns/m
m1=2;      % in kg
m2=4;      % in kg
n=4;        %Number ofStqte vqriqbles

%% System modelling
% State matrix:
A=[0 1 0 0                         
   -(k1+k2)/m1 -alpha/m1 k2/m1 0
   0 0 0 1
   k2/m2 0 -k2/m2 -alpha/m2];

% Control matrix:
B=[0 0 0 1/m2]'; 

% Output matrices:
C1=[1 0 0 0] ;                       
C2=[0 0 1 0];
%Complete C3:
C3=[-1 0 1 0];

% Output/input coupling matrix 
D=zeros(2,1);                      


%% QUESTION 1: Stability, controllability and observability

fprintf("\n=== Q1 Stability ===\n");
eigA = eig(A);
disp("Eigenvalues(A) = "); disp(eigA);

if all(real(eigA) < 0)
    fprintf("System is STABLE (all eigenvalues have negative real part).\n");
else
    fprintf("System is NOT stable.\n");
end

% Transfer functions for cases 1–3 
sys1 = ss(A,B,C1,0);
sys2 = ss(A,B,C2,0);
sys3 = ss(A,B,C3,0);

G1 = tf(sys1); G2 = tf(sys2); G3 = tf(sys3);

fprintf("\n=== Q1: Transfer functions ===\n");
disp("Case 1: G1(s) = X1(s)/U(s)                "); 
disp("Case 2: G2(s) = X3(s)/U(s)                "); 
disp("Case 3: G3(s) = Y_(X3-X1)(s)/U(s)         "); 


fprintf("\n=== Q1: Controllability ===\n");
Rc = rank(ctrb(A,B));
fprintf("rank(ctrb(A,B)) = %d (need 2)\n", Rc);
fprintf("Controllable? %s\n", string(Rc==n));

fprintf("\n=== Q1: Observability ===\n");
Ro1 = rank(obsv(A,C1));
Ro2 = rank(obsv(A,C2));
Ro3 = rank(obsv(A,C3));

fprintf("Case 1 (C1): rank = %d -> Observable? %s\n", Ro1, string(Ro1==n));
fprintf("Case 2 (C2): rank = %d -> Observable? %s\n", Ro2, string(Ro2==n));
fprintf("Case 3 (C3): rank = %d -> Observable? %s\n", Ro3, string(Ro3==n));

%% QUESTION 2: Pole placement - Luenberger observer

pole=[-0.5+i -0.5-i -1 -1.5];

% We must pick which measurement is used by the observer.
pole_sets = {pole};

% First sensor used in the observer 
% Observer as a state-space system (if possible):
C_meas = C1;
L_sets = cell(size(pole_sets));

fprintf("\n=== Q2: Observer gains L for pole set (using C1) ===\n");
for i=1:numel(pole_sets)
    poles = pole_sets{i};
    L_sets{i} = place(A', C_meas', poles)'; % duality
    fprintf("Poles = %s\n", mat2str(poles));
    disp("L = "); disp(L_sets{i});
    disp("eig(A - L*C) = "); disp(eig(A - L_sets{i}*C_meas));
end

% Steady State Representation and transfer functions for case 1
sys1 = ss(A,B,C1,0)
G1 = tf(sys1)



% Second sensor used in the observer
% Observer as a state-space system (if possible):
C_meas = C2;
L_sets = cell(size(pole_sets));

fprintf("\n=== Q2: Observer gains L for pole set (using C2) ===\n");
for i=1:numel(pole_sets)
    poles = pole_sets{i};
    L_sets{i} = place(A', C_meas', poles)'; % duality
    fprintf("Poles = %s\n", mat2str(poles));
    disp("L = "); disp(L_sets{i});
    disp("eig(A - L*C) = "); disp(eig(A - L_sets{i}*C_meas));
end

% Steady State Representation and transfer functions for case 2
sys2 = ss(A,B,C2,0)
G2 = tf(sys2)

% Third sensor used in the observer
% Observer as a state-space system (if possible):

%Nelson told us to cancel this 3rf number configuration





%% QUESTION 3: Diagnosis: study of the transfer functions between the residuals and the faults.

% First sensor used in the observer (j=1)
Cj = C1;
Lj = L_sets{1};   % observer gain computed previously
% State space representation of the state estimation error
Ar = A - Lj*Cj;
Br = -Lj;
Cr = Cj;
Dr = eye(size(Cj,1));

sys_r1 = ss(Ar, Br, Cr, Dr);

fprintf("\n=== Q3: Residual dynamics (Sensor 1) ===\n");
disp("Ar = "); disp(Ar);
disp("Br = "); disp(Br);
disp("Cr = "); disp(Cr);
disp("Dr = "); disp(Dr);

% Transfer functions residual / fault
[num1, den1] = ss2tf(Ar, Br, Cr, Dr, 1);

fprintf("Transfer function r1(s)/f1(s):\n");
tf(num1, den1)

% Second sensor used in the observer (j=2)
Cj = C2;
Lj = L_sets{1};   % observer gain computed previously
% State space representation of the state estimation error
Ar = A - Lj*Cj;
Br = -Lj;
Cr = Cj;
Dr = eye(size(Cj,1));
sys_r2 = ss(Ar, Br, Cr, Dr);


fprintf("\n=== Q3: Residual dynamics (Sensor 2) ===\n");
disp("Ar = "); disp(Ar);
disp("Br = "); disp(Br);
disp("Cr = "); disp(Cr);
disp("Dr = "); disp(Dr);


% Transfer functions residual / fault
[num2, den2] = ss2tf(Ar, Br, Cr, Dr, 1);

fprintf("Transfer function r2(s)/f2(s):\n");
tf(num2, den2)


%Thanks for the luck wishes

