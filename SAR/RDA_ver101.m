% 2016/10/23
%Liu Yakun
clc;clear;close all;
C = 3e8;
Fc = 2e9;
lambda = C / Fc;
H = 3000;
Theta_center = 45 / 180 * pi;
Yc = H * tan(Theta_center);
V = 150;
R0 = sqrt(H^2 + Yc^2);

X_length = 150;
Y_length = 100;
Xmin = 0;
Xmax = X_length;
Ymin = Yc - Y_length/2;
Ymax = Yc + Y_length/2;

La = 2;
Lsar = lambda / La * R0;
Tsar = Lsar / V;
Rmin = sqrt(H^2 + Ymin^2);
Rmax = sqrt(H^2 + Ymax^2 + (Lsar/2)^2);

Tr = 1.5e-6;
Br = 150e6;
Kr = Br / Tr;
alpha_Fsr = 1.2;
Fsr = alpha_Fsr * Br;
Tsr = 1 / Fsr;
Nr = (2*(Rmax-Rmin)/C+Tr)/Tsr;
Nr = 2^nextpow2(Nr);
Tfast = linspace(2*Rmin/C,2*Rmax/C+Tr,Nr);
Tsr = ((2*(Rmax-Rmin)/C+Tr))/Nr;

Ka = -2*V^2/lambda/R0;
Ba = abs(Ka * Tsar);
alpha_Fsa = 1.25;
PRF = alpha_Fsa * Ba;
PRT = 1 / PRF;
Na = (X_length + Lsar) / V / PRT;
Na = 2^nextpow2(Na);
Tslow = linspace((Xmin-Lsar/2)/V,(Xmax+Lsar/2)/V,Na);
PRT = (X_length + Lsar) / V / Na;
PRF = 1/PRT;

Targets = [  0    0    1                  
            -20   0    1 
             50   0    1  
              0  20    1   
             20  20    1];

N = size(Targets,1);
scene_center = [0 Yc];
Targets(:,1:2) = Targets(:,1:2) + ones(N,1) * scene_center; 

echo = zeros(Na,Nr);
for i = 1:N
    delta_x = Tslow * V - Targets(i,1);
    delta_y = Targets(i,2);
    range = sqrt(delta_x.^2+delta_y^2+H^2);
    t = 2 * range / C;
    tau = ones(Na,1) * Tfast - t.' * ones(1,Nr);
    phase = pi*Kr*tau.^2 - pi*4/lambda*range.' * ones(1,Nr);
    rcs = Targets(i,3);
    echo = echo + rcs * exp(1i * phase) .* (tau>0 & tau<Tr) .* ((abs(delta_x) < Lsar/2)' * ones(1,Nr));
end
% figure;
% waterfall(abs(echo));
figure;
imagesc(abs(echo));

t = Tfast - 2*Rmin/C;
refr = exp(1i*Kr*pi*t.^2) .* (t>0 & t<Tr);
signal_comp = ifty(fty(echo) .* (ones(Na,1) * conj(fty(refr))));
figure;
colormap(gray);
imagesc(255-abs(signal_comp));

signal_rD = ftx(signal_comp);
win = waitbar(0,'最近邻域插值');
for i = 1:Na
    for j = 1:Nr
        delta_R = 1/8*(lambda/V)^2*((j-Nr/2)*Tsr*C/2+R0)*((i-Na/2)/Na*PRF)^2;
        nRCM = delta_R * 2 / C / Tsr;
        frac_nRCM = nRCM - floor(nRCM);
        
        if j + round(nRCM) > Nr
            signal_rD(i,j) = signal_rD(i,Nr/2);
        else
            if frac_nRCM < 0.5
                signal_rD(i,j) = signal_rD(i,j+floor(nRCM));
            else
                signal_rD(i,j) = signal_rD(i,j+floor(nRCM)+1);
            end
        end
    end
    waitbar(i/Na);
end
close(win);
% signal = zeros(Na,Nr);
% win = waitbar(0,'最近邻域插值');
% for i = 1:Na
%     for j = 4:Nr
%         delta_R = 1/8*(lambda/V)^2*((j-Nr/2)*Tsr*C/2+R0)*((i-Na/2)/Na*PRF)^2;
%         nRCM = delta_R * 2 / C / Tsr;
%         frac_nRCM = nRCM - floor(nRCM);
%         
%         n_core = 4;
%         for k = -n_core/2:n_core/2-1
%             if j+ceil(nRCM)+k > Nr
%                 signal(i,j) = signal(i,j) + signal_rD(i,Nr)*sinc(k+nRCM);
%             else
%                 signal(i,j) = signal(i,j) + signal_rD(i,j+floor(nRCM)+k)*sinc(k+frac_nRCM);
%             end
%         end
%     end
%     waitbar(i/Na);
% end
% close(win);



signal_rmc = iftx(signal_rD);
figure;
colormap(gray);
imagesc(255-abs(signal_rmc));

t = Tslow - Xmin/V;
refa = exp(1i*pi*Ka*t.^2) .* (abs(t)<Tsar/2);
final_signal = iftx(ftx(signal_rmc) .* (conj(ftx(refa)).' * ones(1,Nr)));
figure;
colormap(gray);
imagesc(255-abs(final_signal));
                
           