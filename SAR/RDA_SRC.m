%2016/12/1
%Liu Yakun

clc,clear,close all;
%%%%%%%%%%%%%%%%%%%%%%参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 3e8;                  %光速
Theta = 45 / 180 * pi;    %斜视角
Rc = 41.7e3;              %场景中心斜距
R0 = Rc * cos(Theta);     %场景中心最近距离
lambda = 0.03;            %发射信号波长
fc = C / lambda;          %发射信号中心频率
Tr = 2e-6;                %脉冲 持续时间
Br = 60e6;                %带宽
Kr = Br / Tr;             %距离调频率
V = 200;                  %平台飞行速度
D = 5;                    %方位向天线长度
Beta = lambda / D;        %波束宽度
Lsar = lambda * Rc / D;   %正侧视合成孔径长度
Lsar_squint = R0 *(tan(Theta + Beta/2) - tan(Theta - Beta/2));
Tsar = Lsar_squint / V;
Dx = D / 2;
Dy = C / 2 / Br;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%场景参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scene_center = [Rc*cos(Theta) Rc*sin(Theta)];
range_width = 500;
azimuth_width = 500;
Ymin = scene_center(1) - range_width/2;
Ymax = scene_center(1) + range_width/2;
Xmin = scene_center(2) - azimuth_width/2;
Xmax = scene_center(2) + azimuth_width/2;


targets = [ scene_center(1)       scene_center(2)      0.8
            scene_center(1)+200   scene_center(2)      0.5
%             scene_center(1)-200   scene_center(2)      0.3
%             scene_center(1)       scene_center(2)+200  1
%             scene_center(1)       scene_center(2)-200  1
%             scene_center(1)+200   scene_center(2)+200  0.7
%             scene_center(1)-200   scene_center(2)-200  1
%             scene_center(1)+200   scene_center(2)-200  0.9
%             scene_center(1)-200   scene_center(2)+200  0.4
            ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%显示场景中的点坐标%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(targets(:,1),targets(:,2),'b*');
xlabel('Range axis [m]');
ylabel('Azimuth axis [m]');
title('Targets in scene');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%慢时间%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ta_start = (Xmin - Ymax * tan(Theta + Beta/2)) / V;
ta_end = (Xmax - Ymin * tan(Theta - Beta/2)) / V;

Ka = -2*V^2*cos(Theta)^2/lambda/Rc;
Ba = abs(Ka*Tsar);
alpha_PRF = 1.3;
PRF = round(alpha_PRF * Ba);
PRT = 1 / PRF;
Na = (ta_end - ta_start) * PRF; 
Na = 2^nextpow2(Na);

fdoc = 2 * V * sin(Theta) / lambda;
% Mamb = floor(fdoc / PRF);
% fdoc = fdoc - Mamb*PRF;
Tslow = [-Na/2:Na/2 - 1] * PRT;
distance_azimuth = Tslow * V;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%快时间%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rnear = Ymin / cos(Theta - Beta/2);
Rfar = Ymax / cos(Theta + Beta/2);

alpha_Fr = 1.2;
Fr = round(alpha_Fr * Br);
Nr = (2*(Rfar - Rnear)/C + Tr) * Fr;
Nr = 2^nextpow2(Nr);

Tfast = [-Nr/2:Nr/2 - 1]/Fr + 2 * Rc / C;
distance_range = Tfast * C / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%回波产生%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo = zeros(Na,Nr);
nTar = size(targets,1);

for i = 1:nTar
    range = sqrt(targets(i,1)^2 + (targets(i,2) - V*Tslow).^2);
    tau = 2 * range / C;
    D = ones(Na,1) * Tfast - tau.' * ones(1,Nr);
    phase = pi*Kr*D.^2 - 4*pi/lambda*(range.' * ones(1,Nr));
    sigma = targets(i,3);
    tci = (targets(i,2) - scene_center(2) + tan(Theta) * (targets(i,1) - scene_center(1))) / V;
    tsi = targets(i,1) * (tan(Theta + Beta/2) - tan(Theta)) / V;
    tei = targets(i,1) * (tan(Theta) - tan(Theta - Beta/2)) / V;
%     time_wave_center = (targets(i,2) - scene_center(2))/V;
    echo = echo + sigma * exp(1i * phase) .* (abs(D) < Tr/2) .* (((Tslow > (tci - tsi) & Tslow < (tci + tei)))' * ones(1,Nr));
end

figure;
imagesc(abs(echo));
xlabel('Range axis');
ylabel('Azimuth axis');
title('Raw signal');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%Range Chirp signal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% chirp_range = zeros(1,Nr);
% t_range = -Tr/2:1/Fr:Tr/2;
% omega_range = -Fr/2:1/Tr:Fr/2;
% chirp_range_temp = exp(1i * pi * Kr * t_range.^2);
% length_range = length(t_range);
% chirp_range(ceil((Nr - length_range) / 2):ceil((Nr - length_range) / 2) + length_range - 1) = chirp_range_temp;
% 
% figure;
% subplot(211);
% plot(t_range,real(chirp_range_temp),'r');
% hold on;
% plot(t_range,imag(chirp_range_temp),'g');
% xlabel('Time axis [sec]');
% ylabel('Amplitude');
% xlim([min(t_range) max(t_range)]);
% title('Range Chirp');
% 
% subplot(212);
% plot(omega_range,abs(fftshift(fft(chirp_range_temp))));
% xlabel('Frequence [Hz]');
% ylabel('Absolute');
% xlim([min(omega_range) max(omega_range)]);
% title('Spectrum of range chirp');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Azimuth Chirp signal%%%%%%%%%%%%%%%%%%%%%%
% chirp_azimuth = zeros(1,Na);
% t_azimuth = -Tsar/2:1/PRF:Tsar/2;
% omega_azimuth = (fdoc - PRF/2):1/Tsar:(fdoc + PRF/2);
% chirp_azimuth_temp = exp(1i * pi * Ka * t_azimuth.^2);
% length_azimuth = length(t_azimuth);
% chirp_azimuth(ceil((Nr - length_azimuth) / 2):ceil((Nr - length_azimuth) / 2) + length_azimuth - 1) = chirp_azimuth_temp;
% 
% figure;
% subplot(211);
% plot(t_azimuth,real(chirp_azimuth_temp),'r');
% hold on;
% plot(t_azimuth,imag(chirp_azimuth_temp),'g');
% xlabel('Time axis [sec]');
% ylabel('Amplitude');
% xlim([min(t_azimuth) max(t_azimuth)]);
% title('Azimuth Chirp');
% 
% subplot(212);
% plot(omega_azimuth,abs(fftshift(fft(chirp_azimuth_temp))));
% xlabel('Frequence [Hz]');
% ylabel('Absolute');
% xlim([min(omega_azimuth) max(omega_azimuth)]);
% title('Spectrum of range chirp');
%%%%%%%%%%%%%%%%%%%%%%距离压缩%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = Tfast - 2 * Rc / C;
chirp_range = exp(1i * pi * Kr * t.^2) .* (abs(t) < Tr / 2);
signal_Ra = fty(echo) .* (ones(Na,1) * conj(fty(chirp_range))); 
signal_comp = ifty(signal_Ra);

%%%%%%%%%%%%%%%%%%%%%%%%%%显示距离压缩后的图像%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
colormap('gray');
imagesc(255 - abs(signal_comp));
xlabel('Range axis');
ylabel('Azimuth axis');
title('Range compression');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SRC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal_RA = ftx(signal_Ra);
% f_azimuth = zeros(1,Na);
% f_azimuth_temp = (fdoc - PRF/2):(1 / Total_time):(fdoc + PRF/2);
% length_azimuth = length(t_azimuth);
% % f_azimuth(ceil((Na - length_azimuth)/2):ceil((Na - length_azimuth)/2) + length_azimuth - 1) = f_azimuth_temp;
% f_range = zeros(1,Nr);
% f_range_temp= -Fr/2:1/Tr:Fr/2;
% length_range = length(t_range);
% f_range(ceil((Nr - length_range) / 2):ceil((Nr - length_range) / 2) + length_range - 1) = f_range_temp;
f_range = linspace(-Fr/2,Fr/2,Nr);
f_azimuth = linspace(fdoc - PRF/2,fdoc + PRF/2,Na);
D = sqrt(1 - (lambda * f_azimuth / 2 / V).^2);
Ksrc = 2 * V^2 * fc^3 / C / R0 .* D.^3 ./f_azimuth;

H_src = exp(-1i * pi * (ones(Na,1) * f_range.^2) ./(Ksrc' * ones(1,Nr)));
signal_RA_src = signal_RA .* H_src;
signal_rA_src = ifty(signal_RA_src);
signal_ra_src = iftx(signal_rA_src);
signal_ra = iftx(ifty(signal_RA));

figure;
subplot(211);
% colormap('gray');
imagesc(abs(signal_ra));
title('No SRC');

subplot(212);
% colormap('gray');
imagesc(abs(signal_ra_src));
title('After SRC');
signal_RD = signal_rA_src;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RCMC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
win = waitbar(0,'Sinc 插值');
length_core = 4;
signal_rcmc = zeros(Na,Nr);

for i = 1:Na
    for j = 1:Nr
        fa = f_azimuth(i);
        D = sqrt(1 - (fa*lambda/2/V)^2);
        r = (Rc + (i - Nr/2)/Fr * C / 2) * cos(Theta);
        rcm = r * (1/D - 1);
        n_rcm = 2 * r / C * Fr;
        delta_nrcm = n_rcm - floor(n_rcm);
        
        if j+n_rcm > Nr
            signal_rcmc(i,j) = signal_RD(i,Nr/2);
        else
            if delta_nrcm < 0.5
                signal_rcmc(i,j) = signal_RD(i,j + floor(n_rcm));
            else
                signal_rcmc(i,j) = signal_RD(i,j + ceil(n_rcm));
            end
        end
    
%         for k = -length_core/2 : length_core/2 - 1
%             if j+rcm+k > Nr
%                 signal_rcmc = signal_rcmc + signal_RD(i,Nr) * sinc(k+n_rcm);
%             else
%                 signal_rcmc = signal_rcmc + signal_RD(i,j + round(n_rcm) + k) * sinc(delta_nrcm + k);
%             end
%         end
    end
    waitbar(i/Na);
end
   
close(win);

figure;
imagesc(abs(iftx(signal_rcmc)));
title('After RCMC');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Azimuth Compression%%%%%%%%%%%%%%%%%%%%%%%%%
t = Tslow - scene_center(2)/V;
chirp_azimuth = exp(1i * pi * Ka * t.^2) .* (abs(t) < Tsar/2);
signal_RD = ftx(signal_comp);
signal_final = iftx(signal_RD .* (conj(ftx(chirp_azimuth)).' * ones(1,Nr)));

figure;
imagesc(abs(signal_final));
title('Final Signal');