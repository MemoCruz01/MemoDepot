clear all
close all
clc;
p=tf('s')
%% h=0 et syst�me sans retard

% close all
G=10/(p*(1+0.16*p));
G.IODelay=0;

w=logspace(-4,4,100);

% Etude 
figure;
nichols(G,w)
ngrid
hold all
%% Question 1: Gain �quivalent de la non lin�arit�
% La non-lin�arit� s'�crit N(x1,w)=rho*exp(-jalpha)
%avec rho=4*M/(pi*x1)
% et alpha=arcsin(h/(2*x1))
%% Question 2
%% Partie 1
condition_initiale_init=0;
h=0;
retard=0;
[t,x,y] = sim('Exo_NL_1_Activite_finale')
consigne=y(:,1);
mesure=y(:,2);
erreur=y(:,3);
erreur_sat=y(:,4);

figure;
ax(1)=subplot(211)
plot(t,mesure,t,consigne)
legend('mesure','consigne')
ax(2)=subplot(212)
plot(t,erreur_sat,t,erreur)
legend('Erreur satur�e','Erreur')

%% Partie 2: Trac� des lieux de transfert dans le plan de Nichols
x1=[0.00001:0.001:10000];
M=pi;
N1=4*M./(pi*x1);

figure;
plot(x1,N1,'o')
legend('N(x1)')
figure;
plot(x1,-1./N1,'x')
legend('-1/N(x1)')

% Calcul de gain et de phase de la NL et de la B.O.
Gain_NL=20*log10(1./N1);
figure;
plot(-180.*[ones(length(Gain_NL),1)],Gain_NL,'linewidth',3)
[mag,phase]=bode(G,w)
hold all
mag_db(1:length(mag))=20*log10(mag(1,1,1:length(mag)));
phase_degre(1:length(phase))=phase(1,1,1:length(phase));
plot(phase_degre,mag_db)
ngrid
%==> Aucune intersection entre les lieux de transfert donc il ne devrait
%pas apparaitre d'oscillations (confirm� par le trac� temporelle de la
%% Question 3
%% partie 1 : h=0.06 et syst�me sans retard

G=10/(p*(1+0.16*p));
G.IODelay=0;

%Partie 1
condition_initiale_init=0;
h=0.06;
retard=0;
M=pi;
[t,x,y] = sim('Exo_NL_1_Activite_finale')
consigne=y(:,1);
mesure=y(:,2);
erreur=y(:,3);
erreur_sat=y(:,4);
%==> apparait une oscillation d'amplitude (0.17*2)=0.034 et de fr�quence 6Hz
figure;
ax(1)=subplot(211)
plot(t,mesure,t,consigne)
legend('mesure','consigne')
ax(2)=subplot(212)
plot(t,erreur_sat,t,erreur)
legend('Erreur satur�e','Erreur')

%% Partie 2: Trac� des lieux de transfert dans le plan de Nichols
% Rappel : La non-lin�arit� s'�crit N(x1,w)=rho*exp(-jalpha)
%avec rho=4*M/(pi*x1)
% et alpha=arcsin(h/(2*x1))

w=logspace(-4,4,100)
x1=[(h)/2:0.001:10000];
N1=(4*M./(pi*x1)).*exp(-i.*asin(h./(2*x1)));
figure;
plot(x1,N1,'o')
legend('N(x1)')
figure;
plot(x1,-1./N1,'x')
legend('-1/N(x1)')

%% Calcul de gain et de phase du NL
Gain_NL=20*log10(1./(4*M./(pi*x1)));
figure;
nichols(G)
hold all
% G1=(10*0.1)/(p*(1+0.16*p));
% G1.IODelay=0;
% nichols (G1)
plot(-180+180*asin(h./(2*x1))/pi,Gain_NL,'linewidth',3);%phase= -180+dephasage de NL
% [mag,phase]=bode(G,w)
% hold all
% mag_db(1:length(mag))=20*log10(mag(1,1,1:length(mag)));
% phase_degre(1:length(phase))=phase(1,1,1:length(phase));
% plot(phase_degre,mag_db)
ngrid
% Point d'intersection a:
gain_db=-26.94;
x10=4*10^(gain_db/20)
f10=37.1/(2*pi)

%% Question 4 : h=0 et syst�me a retard 

close all
G=10/(p*(1+0.16*p));
G.IODelay=0.02;
M=pi;
w=logspace(-4,4,100)

% Etude 
figure;
nichols(G,w)
ngrid
hold all

% Gain non lin�arit�
x1=[0.00001:0.00001:1];
N1=4*M./(pi*x1);
figure;
plot(x1,N1,'o')
figure;
plot(x1,-1./N1,'x')

%% Partie 1
condition_initiale_init=0;
h=0;
retard=0.02;
[t,x,y] = sim('Exo_NL_1_Activite_finale')
consigne=y(:,1);
mesure=y(:,2);
erreur=y(:,3);
erreur_sat=y(:,4);
%==> apparait une oscillation 
figure;
ax(1)=subplot(211)
plot(t,mesure,t,consigne)
legend('mesure','consigne')
ax(2)=subplot(212)
plot(t,erreur_sat,t,erreur)
legend('Erreur satur�e','Erreur')

%% Partie 2: Trac� des lieux de transfert dans le plan de Nichols
% Rappel : La non-lin�arit� s'�crit N(x1,w)=rho*exp(-jalpha)
%avec rho=4*M/(pi*x1)
% et alpha=arcsin(h/(2*x1))

Gain_NL=20*log10(1./N1);
figure;
plot(-180.*[ones(length(Gain_NL),1)],Gain_NL,'linewidth',3);
[mag,phase]=bode(G,w);
hold all
mag_db(1:length(mag))=20*log10(mag(1,1,1:length(mag)));
phase_degre(1:length(phase))=phase(1,1,1:length(phase));
plot(phase_degre,mag_db);
hold all % Altrenative au trac� de nichols de G
nichols(G)
ngrid

%% Calcul de x10 
gain_db=-14.15;
x10=4*10^(gain_db/20)
f10=17.4/(2*pi)
%% Question 5
% Calcul de w=1/tau
close all
G=10/(p*(1+0.16*p));
G.IODelay=0.02;
M=pi;
% Exemple : Avance de phase
a=5
tau=10
AP=(1+(a*tau*p))/(1+(tau*p/a))
figure;
bode(AP)
%%
% L'avance phase apport�e par C(p) est maximale � 1/tau et vaut:
a=5; % Donn� par l'�nonc�
PhiMax=2*180*atan(a)/pi-90
% On souhaite que la pulsation de pompage w0 soit la plus grande possible donc
% on choisira w0=1/Tau
%==>�quation permettant de choisir w0 : arg(G(jw0)+PhiMax=-180
%==> on recherche donc w0 tel que arg(G(jw0)=-180-PhiMax
Arg_rech=-180-PhiMax
Gain_NL=linspace(-100,100,100);
figure;
plot(Arg_rech.*[ones(length(Gain_NL),1)],Gain_NL,'linewidth',3);
hold all
nichols(G)
ngrid
w0=63.7;
Tau=1/w0
C=(1+5*Tau*p)/(1+Tau*p/5);
figure;
Gain_NL=20*log10(1./N1);
plot(-180.*[ones(length(Gain_NL),1)],Gain_NL,'linewidth',3);
hold all
nichols(G)
hold all
nichols(G*C)