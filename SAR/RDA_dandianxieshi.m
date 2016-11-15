%2016/11/11
%Liu Yakun

clc;
clear;
close all;

C = 3e8;
Fc = 3e9;
lambda = C / Fc;
V = 150;
Theta = 10 / 180 * pi;
Yc = 10000;
D = 4;
beta = lambda / D;
Xc = Yc * tan(Theta);
Rc = Yc / cos(Theta);

Tslow_begin = -Yc * (tan(Theta + beta/2) - tan(Theta)) / V;
Tslow_end = Yc * (tan(Theta) - tan(Theta - beta/2)) / V;
Tsar = Tslow_end - Tslow_begin;
Ka = 2 * V^2 * cos(Theta)^2 / lambda / Rc;
Ba = Ka * Tsar;
alpha_PRF = 1.3;
alpha_Fs = 1.2;
PRF = round(alpha_PRF * Ba);
Na = (Tslow_end - Tslow_begin) / V * PRF;

