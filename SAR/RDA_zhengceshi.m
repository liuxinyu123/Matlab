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
PRT = 1 / PRF;
Na = round((Ta_end - Ta_begin) * PRF);
Na = 2^nextpow2(Na);
Tslow = [-Na/2:Na/2 - 1] * PRT;

Fr_alpha = 1.2;
Fs = Fr_alpha * Br;
Ts = 1 / Fs;
Rnear = Yc - Wr/2;
Rfar = (Yc + Wr/2)/cos(beta/2);
Rmid = (Rnear + Rfar) / 2;
Nr = (2 * (Rfar - Rnear) / C + Tr)/Ts;
Nr = 2^nextpow2(Nr);

Tfast = [-Nr/2:Nr/2-1]*Ts + (2 * Rmid / C);

ntarget = size(Targets,1);
echo = zeros(Na,Nr);

for i = 1:ntarget
    x = Targets(i,1);
    y = Targets(i,2);
    rcs = Targets(i,3);
    range = sqrt(y^2 + (x - V*Tslow).^2);
    tau = 2 * range / C;
    t_center = x / V;
    delta_t = y * tan(beta/2);
    D = ones(Na,1)*Tfast - tau' * ones(1,Nr);
    phase = pi * D.^2 - 4 * pi / lambda * range' * ones(1,Nr);
    t = ones(Na,1) * Tfast - tau' * ones(1,Nr);
    echo = echo + rcs*exp(1i*phase) .* ((abs(Tslow - t_center) < delta_t)' * ones(1,Nr)) .* ((abs(t) < Tr/2));
end

t = Tfast - 2 * Rmid / C;
ref_r = exp(1i*pi*Kr*t.^2) .* (abs(t) < Tr/2);% .* hamming(Nr).';
ref_R = fty(ref_r);
signal_comp = ifty(fty(echo) .* (ones(Na,1) * conj(ref_R)));


