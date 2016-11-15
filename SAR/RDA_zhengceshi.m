%2016/11/14
%Liu Yakun

clc,clear,close all;

C = 3e8;
Fc = 1e9;
lambda = C / Fc;
V = 150;
D = 2;
beta = lambda / D;
Br = 60e6;
Tr = 2e-6;
Kr = Br / Tr;
Yc = 10e3;
Xc = 0;
Lsar = beta * Yc;
Tsar = Lsar / V;
Ka = 2 * V^2 / lambda / Yc;
Ba = abs(Ka * Tsar);
Wr = 300;
Wa = 200;
a = 0.3;
b = 0.4;
Targets = [Xc + a * Wa/2,Yc + b * Wr/2,1
           Xc + a * Wa/2,Yc - b * Wr/2,1
           Xc - a * Wa/2,Yc + b * Wr/2,1
           Xc - a * Wa/2,Yc - b * Wr/2,1
           Xc           ,Yc           ,1];
       
Ta_begin = -((Yc + Wr/2)*tan(beta/2) + Wa/2) / V;
Ta_end = ((Yc + Wr/2)*tan(beta/2) + Wa/2) / V;
Fa_alpha = 1.3;
PRF = Fa_alpha * Ba;
Na = round((Ta_end - Ta_begin) * PRF);
Na = 2^nextpow2(Na);
PRT = (Ta_end - Ta_begin) / Na;
Tslow = [-Na/2:Na/2 - 1] * PRT;


