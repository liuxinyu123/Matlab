%2016/12/15
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
range_width = 200;
azimuth_width = 500;
Ymin = scene_center(1) - range_width/2;
Ymax = scene_center(1) + range_width/2;
Xmin = scene_center(2) - azimuth_width/2;
Xmax = scene_center(2) + azimuth_width/2;


targets = [ scene_center(1)       scene_center(2)        1
%             scene_center(1)+100   scene_center(2)        1
%             scene_center(1)-100   scene_center(2)        1
            scene_center(1)       scene_center(2)+100  1
            scene_center(1)       scene_center(2)-100  1
%             scene_center(1)+100   scene_center(2)+100  0.7
%             scene_center(1)-100   scene_center(2)-100  1
%             scene_center(1)+100   scene_center(2)-100  0.9
%             scene_center(1)-100   scene_center(2)+100  0.4
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
    phase = -pi*Kr*D.^2 - 4*pi/lambda*(range.' * ones(1,Nr));
    sigma = targets(i,3);
    tci = (targets(i,2) - scene_center(2) + tan(Theta) * (targets(i,1) - scene_center(1))) / V;
    tsi = targets(i,1) * (tan(Theta + Beta/2) - tan(Theta)) / V;
    tei = targets(i,1) * (tan(Theta) - tan(Theta - Beta/2)) / V;
%     time_wave_center = (targets(i,2) - scene_center(2))/V;
    echo = echo + sigma * exp(1i * phase) .* (abs(D) < Tr/2) .* (((Tslow > (tci - tsi) & Tslow < (tci + tei)))' * ones(1,Nr));
end

% R_ref = Rc * cos(Theta);%参考距离
% fa = [-Na/2:Na/2 - 1]/Na * PRF;%方位向频率
% fr = [-Nr/2:Nr/2 - 1]/Nr * Fr;%距离向频率
% D = sqrt(1-(lambda*fa/2/V).^2);%徙动系数
% D_ref = sqrt(1 - (lambda*fdoc/2/V)^2);%参考方位频率处的徙动系数
% tau = ones(Na,1) * Tfast - 2*R_ref/C./(D' * ones(1,Nr));%距离向时间转换
% R = R_range * cos(Theta);
% Ksrc = 2*V^2*fc^3/C*(D'.^3 ./ fa'.^2 * ones(1,Nr)) ./ (ones(Na,1) * R_range);
% Km = Kr ./ (1 - Kr./Ksrc);
% Ksrc_ref = 2*V^2*fc^3/C/R_ref * (D.^3./fa.^2);
% Ks = Kr ./ (1 - Kr./Ksrc_ref);
% 
% signal_rA = ftx(echo);
% signal_RA = fty(signal_rA);
% three_phase = pi*C/2/V^2/fc^4*(fa'.^2./D'.^5 * ones(1,Nr)) .* (ones(Na,1) *  (R .* fr.^3));
% S_three_cor = exp(1i * three_phase);
% signal_RA_cor = signal_RA .* S_three_cor;
% signal_rA_cor = ifty(signal_RA_cor);
% 
% Ssc_phase = pi * (Ks' * ones(1,Nr)) .* ((D_ref./D - 1)' * ones(1,Nr)) .* tau.^2;
% Ssc = exp(1i * Ssc_phase);
% signal_rA_sc = signal_rA_cor .* Ssc;
% 
% signal_RA_sc = fty(signal_rA_sc);
% range_src_phase = pi * (D'/D_ref./Ks' * ones(1,Nr)) * (ones(Na,1) * fr.^2);
% bulk_rcm_phase = 4*pi/C*((D_ref./D' - 1) * ones(1,Nr)) .* (ones(Na,1) * (R_ref * fr));
% S_range = exp(1i * range_src_phase);
% S_bulk = exp(1i * bulk_rcm_phase);
% signal_RA_range = signal_RA_sc .* S_range;
% signal_RA_bulk = signal_RA_range.* S_bulk;
% signal_RD = ifty(signal_RA_bulk);
% 
% azimuth_phase = 4*pi/lambda*(ones(Na,1) * R) .* (D' * ones(1,Nr));
% S_azimuth = exp(1i * azimuth_phase);
% signal_RD_azimuth = signal_RD .* S_azimuth;
f_doc = 2*V*sin(Theta)/lambda;             %多普勒中心频率
f_range = [-Nr/2:Nr/2-1]/Nr * Fr + fdoc;     %距离向频率
f_azimuth = [-Na/2:Na/2-1]/Na * PRF;   %方位向频率
R = C / 2 * Tfast;                    % 斜距
r = R * cos(Theta);                   %最近斜距
R_ref = Rc;                           %参考距离
D = sqrt(1 - (lambda * f_azimuth/2/V).^2);    %徙动因子
signal_rA = ftx(echo);     %方位向傅立叶变换
signal_RA = fty(signal_rA);  %距离向傅立叶变换
%三次相位补偿
H_ThreePhase = exp(1i*pi*C/2/V^2/fc^4 * ((f_azimuth.^2./D.^5)' * ones(1,Nr)) .* (ones(Na,1)*(r.*f_range))); 
signal_RA = signal_RA .* H_ThreePhase;

signal_rA = ifty(signal_RA);        %距离向傅立叶反变换
Cs_f = cos(Theta)./D - 1;           %CS因子
tau = ones(Na,1) * Tfast -  ((2*R_ref/C*(1+Cs_f))' * ones(1,Nr));         %新的时间
Ksrc = 2*V^2*fc^3/C/(R_ref*cos(Theta)) * ((D.^3./f_azimuth.^2));  %SRC距离调频率
Km = Kr ./ (1 - Kr./Ksrc);        %新的调频率
H_CS = exp(-1i*pi .* ((Km .* Cs_f)' * ones(1,Nr)) .* tau.^2);  %CS信号，去除补余RCM
signal_rA_cs = signal_rA .* H_CS;
signal_RA_cs = fty(signal_rA_cs);          %距离向傅立叶变换

H_rangecomp = exp(-1i*pi ./((Km .* (1+Cs_f))' * ones(1,Nr)) .* (ones(Na,1) * f_range.^2));   %用于距离压缩
H_bulkrcm = exp(1i*4*pi/C*R_ref * (ones(Na,1) * f_range) .* (Cs_f' * ones(1,Nr)));   %用于校正一致RCM
signal_RA_rangecomp_rcm = signal_RA_cs .* H_rangecomp .* H_bulkrcm;
signal_rA_rangecomp_rcm = ifty(signal_RA_rangecomp_rcm);

H_azimuthcomp = exp(1i*4*pi/lambda* (ones(Na,1) * R * cos(Theta)) .* (D' * ones(1,Nr)));
theta1 = -4*pi/C^2* ((Km .* (1+Cs_f).*Cs_f)' * ones(1,Nr)) .* (ones(Na,1) * (R - R_ref).^2);
theta2 = -2*pi*sin(Theta)/V * (f_azimuth' * ones(1,Nr)) .* (ones(Na,1) * R);
H_otherphase = exp(-1i*(theta1+theta2));
signal_rA_azimuthcomp = signal_rA_rangecomp_rcm .* H_azimuthcomp .* H_otherphase;

signal_final = iftx(signal_rA_azimuthcomp);

figure;
imagesc(abs(signal_final));
title('final signal');

