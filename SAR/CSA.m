%2016/12/13
%Liu Yakun

clc,clear,close all;

%%%%%%%%%%%%%%%%%%%%%%参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 3e8;                  %光速
Theta = 0 / 180 * pi;    %斜视角
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


targets = [ scene_center(1)       scene_center(2)        1
%             scene_center(1)+200   scene_center(2)        1
%             scene_center(1)-200   scene_center(2)        1
            scene_center(1)       scene_center(2)+200  1
            scene_center(1)       scene_center(2)-200  1
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

PRF =  2 * round(alpha_PRF * Ba);
Na = (ta_end - ta_start) * PRF; 
Na = 2^nextpow2(Na);

fdoc = 2 * V * sin(Theta) / lambda;
% Mamb = floor(fdoc / PRF);
% fdoc = fdoc - Mamb*PRF;
Tslow = [-Na/2:Na/2 - 1] / PRF;
distance_azimuth = Tslow * V;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%快时间%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rnear = Ymin / cos(Theta - Beta/2);
Rfar = Ymax / cos(Theta + Beta/2);

PRFmin = Ba + 2*V*Kr*Tr*sin(Theta)/C;
RFmax = 1/(2*Tr+2*(Rfar-Rnear)/C);
Rmid = (Rnear + Rfar) / 2;
alpha_Fr = 1.2;
Fr = round(alpha_Fr * Br);
Nr = (2*(Rfar - Rnear)/C + Tr) * Fr;
Nr = 2^nextpow2(Nr);

Tfast = [-Nr/2:Nr/2 - 1]/Fr + 2 * Rmid / C;
R_range = Tfast * C / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%回波产生%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo = zeros(Na,Nr);
nTar = size(targets,1);
disp('*****************************************');
disp('参数打印:');
disp(['慢时间采样点数为:',poly2str(Na),'点']);
disp(['脉冲重复发射频率为:',poly2str(PRF),'Hz']);
disp(['快时间采样点数为:',poly2str(Nr),'点']);
disp(['快时间采样频率为:',poly2str(Fr),'Hz']);
disp(['方位向合成孔径长度为:',poly2str(Lsar),'m']);
disp(['方位向分辨率为:',poly2str(Dx),'m']);
disp(['距离向分辨率为:',poly2str(Dy),'m']);
disp('*****************************************');
disp('产生回波信号.....');
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

R_ref = Rc;%参考距离
fa = [-Na/2:Na/2 - 1]/Na * PRF;%方位向频率
fr = [-Nr/2:Nr/2 - 1]/Nr * Fr;%距离向频率
D = sqrt(1-(lambda*fa/2/V).^2);%徙动系数
D_ref = sqrt(1 - (lambda*fdoc/2/V)^2);%参考方位频率处的徙动系数
tau = ones(Na,1) * Tfast - 2*R_ref/C./(D' * ones(1,Nr));%距离向时间转换
Ksrc = 2*V^2*fc^3/C*(D'.^3 ./ fa'.^2 * ones(1,Nr)) ./ (ones(Na,1) * R_range);
Km = Kr ./ (1 - Kr./Ksrc);
Ssc = exp(1i*pi*Km .* ((D_ref./D - 1)' * ones(1,Nr)) .* tau.^2);%补余RCM参数

signal_rA = fftshift(fft(fftshift(echo))) .* Ssc;%变换到距离多普勒域与补余RCM相位相乘
signal_RA = fftshift(fft(fftshift(signal_rA).')).';%变换到二维频域

range_match_phase = pi/D_ref *(D' * ones(1,Nr)) .* (ones(Na,1) * fr.^2) ./ Km; %距离压缩相位
S_range = exp(1i*range_match_phase);
signal_RA_comp = signal_RA .* S_range;%二维频域中进行距离压缩
bulk_rcm_phase = 4*pi/C*R_ref * ((1./D - 1/D_ref)' * ones(1,Nr)) .* (ones(Na,1) * fr);%一致RCM相位
S_bulk = exp(1i*bulk_rcm_phase);
signal_RA_bulk = signal_RA_comp .* S_bulk; %一致RCM校正

signal_RD = fftshift(ifft(fftshift(signal_RA_bulk).')).';%变换到距离多普勒域
signal_RD1 = fftshift(ifft(fftshift(signal_RA_comp).')).';%对比

azimuth_match_phase = 4*pi*fc/C*(ones(Na,1) * R_range) .* (D' * ones(1,Nr));%方位压缩相位
S_azimuth = exp(1i*azimuth_match_phase);
signal_RD_comp = signal_RD .* S_azimuth;%进行方位压缩
signal_RD_comp1 = signal_RD1 .* S_azimuth;

cor_phase = -4*pi/C^2*Km .* ((1 - D'/D_ref) * ones(1,Nr)) .* (ones(Na,1) * (R_range - R_ref).^2) ./ (D'.^2 * ones(1,Nr));%相位校正相位
S_cor = exp(1i*cor_phase);
signal_RD_cor = signal_RD_comp .* S_cor;%相位校正
signal_RD_cor1 = signal_RD_comp1 .* S_cor;

signal_final = fftshift(ifft(fftshift(signal_RD_cor)));
signal_final1 = fftshift(ifft(fftshift(signal_RD_cor1)));

figure;
imagesc(abs(signal_final));
title('Final Signal');

figure;
plot(abs(signal_final(512,:)));
title('距离向切面');

figure;
plot(10*log(abs(signal_final(512,:)/max(signal_final(512,:)))));
title('峰值旁瓣比');


