clear all 
close all 
clc

s = tf('s');
s_hat = 0.12;
zeta = abs(log(s_hat))/sqrt(pi^2+(log(s_hat))^2);

Ga = 0.014;
Gf = 1;
Gs = 1;
Gp = 100/(s^2+5.5*s+4.5);
kd = 1;
kc = 21.43;
nu = 1;

Lin = (kc)/(s^nu)*Ga*Gf*Gs*Gp

figure(1)
bode(Lin)

[numLin,denLin]=tfdata(Lin,'v');
figure(2)
nyquist1(numLin,denLin),grid on

Tp = 1.078;
Sp = 1.3939;
figure(3)
myngridst(Tp,Sp)
nichols(Lin)

wc= 1.60;
%zero
wnorm = 0.55;
z = wc/wnorm;
zd= 1+s/z;
L = Lin*zd;
figure(3)
nichols(L,'r');

%lead
wnormd = 0.6;
md = 16;
rd = wc/wnormd;
Rd = (1+s/rd)/(1+s/(md*rd));
L = L*Rd;
figure(3)
nichols(L,'g');
%lag
wnorm = 1000;
mi = 10^(8.93/20);
pd = wc/wnorm;
Pd = (1+s/(mi*pd))/(1+(s/pd));
L = L*Pd;
figure(3)
nichols(L,'r');

L = (L*1/(1+(s/60)))*1/(1+(s/60))*1/(1+(s/60));
T = zpk(minreal(L/(1+L),10^-3));
figure(4)
step(T*kd)

figure(5)
bodemag(T)
S = zpk(minreal(1/(1+L),10^-3));
figure,bodemag(S)
