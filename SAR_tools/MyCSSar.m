%=========================================================
% Author: Liu Yakun
% Date: 29,Sep,2016
% Program Name: The image of dot's targets in SAR with CSA
%=========================================================

clc,clear,close all;

%＝＝＝＝＝＝＝＝＝＝＝＝＝参数设置＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
C = 3e+8;       %光速
fc = 5.3e9;     %发射信号中心频率
lambda = C / fc; %波长
Theta_r = 0 / 180 * pi; %斜视角
H = 5000;                  %飞行平台高度
Platform_center = [0,0,H]; %平台坐标
Theta_d = 45 / 180 * pi;   %下视角
Yc = H * tan(Theta_d);     %场景（照射）中心Y轴
Scene_center = [0,Yc,0];   %场景中心坐标
R0 = sqrt(sum((Platform_center - Scene_center).^2));   %天线平台到场景中心的距离
delta_X = 200;
delta_Y = 500;
X_min = -delta_X / 2;
X_max = delta_X / 2;
Y_min = -delta_Y / 2;
Y_max = delta_Y / 2;

%＝＝＝＝＝＝＝＝＝＝＝＝LFM 信号参数＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
Tr = 2.5e-6;       %脉冲持续时间
Kr = 2e+13;        %距离调频率
Br = Tr * Kr;      %发射信号带宽

%雷达性能参数
V = 150;                   %飞行速度
La = 4;                    %天线方位向长度
DY = C / Br / 2;           %距离向分辨率
DX = La / 2;               %方位向分辨率
Lsar = lambda * R0 / La;   %合成孔径长度
Tsar = Lsar / V;           %合成孔径时间

%＝＝＝＝＝＝＝＝＝＝＝快时间参数 距离向＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝

Fr_rate = 1.2;     %距离向过采样系数
Fs = round(Br * Fr_rate); %距离向采样频率
Ts = 1 / Fs;
Rmin = sqrt(H^2 + (Yc + Y_min)^2);                 %最近斜距
Rmax = sqrt(H^2 + (Yc + Y_max)^2 + (Lsar / 2)^2);  %最远斜距
Nr = ceil((2 * (Rmax - Rmin) / C + Tr) / Ts);      %快时间采样点数
Nr = 2 ^ nextpow2(Nr);                             %扩展到2的幂 方便fft运算
Ts = (2 * (Rmax - Rmin) / C + Tr) / Nr;
Fs = ceil(1 / Ts);                                 %更新Fs
tf_ideal = [-Nr / 2 : Nr / 2 - 1] * Ts;              %理想快时间采样时间
% tf = tf_ideal + 2 * R0 / C;                        %实际快时间采样时间
tf = linspace(2 * Rmin / C,2 * Rmax / C + Tr,Nr);  %实际快时间采样时间
tr = tf * C / 2;                                       %快时间采样对应的距离

%＝＝＝＝＝＝＝＝＝＝＝慢时间参数 方位向＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
PRF_rate = 1.25;            %慢时间过采样系数
Ka = -2 * V^2 / lambda / R0; %方位调频率
Ba = abs(Ka * Tsar);              %多普勒带宽
PRF = round(PRF_rate * Ba);  %脉冲重复发射频率
PRT = 1 / PRF;               %脉冲重复发射周期
Na = ceil((X_max - X_min + Lsar) / V / PRT);  %慢时间采样点数
Na = 2 ^ nextpow2(Na);                        %扩展到2的幂 方便fft运算
PRT = (X_max - X_min + Lsar) / V / Na;
PRF = ceil(1 / PRT);                           %更新PRF
ts = [-Na / 2 : Na / 2 - 1] * PRT;                  %慢时间采样时间
ta = ts * V;                                  %慢时间采样对应的距离

%=================目标参数设置＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
Targets = [ 0   0   0   1
           -50 -200 0   1
            50 -200 0   1
            20  200 0   1
           -20  200 0   1];    %［x,y,z,rcs]  相对坐标
Ntar = size(Targets,1);
Targets(:,1:3) = Targets(:,1:3) + ones(Ntar,1) * Scene_center; %实际坐标


%====================回波产生＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
Sr = zeros(Na,Nr);
for i = 1:1:size(Targets,1)
    
    sigma = Targets(i,4);                                  %当前目标的RCS
    x_range = ts * V - Targets(i,1);
    y_range = Platform_center(2) - Targets(i,2);   
    z_range = Platform_center(3) - Targets(i,3);
    range = sqrt(x_range .^ 2 + y_range ^ 2 + z_range ^ 2); %慢时间对应的瞬时斜距
    tau = 2 * range / C;                                    %当前斜距对应的时间延迟
    Dfast = ones(Na,1) * tf - tau' * ones(1,Nr);            %时间矩阵
    phase = pi * Kr * Dfast .^ 2 + 4 * pi / lambda * range' * ones(1,Nr);  %接收信号相位
    Sr = Sr + exp(1i * phase) .* (abs(Dfast) < Tr / 2) .* ((abs(x_range) < Lsar / 2)' * ones(1,Nr)); %接收信号
    
end

%====================CS参数==============================================
fr = linspace(-Fs/2,Fs/2,Nr);       %距离向频率
fa = linspace(-PRF/2,PRF/2,Na);     %方位向频率
R_ref = R0;                         %场景中心最短距离
R = tr;                             %距离向斜距
D = 1 ./ sqrt(1 - (lambda * fa ./ 2 / V).^ 2);
alpha = 1 ./ D - 1;
Km = Kr ./ (1 - Kr * lambda * R_ref * (fa .^ 2) / 2 / V^2 / fc^2 ./ D.^3);
tau = ones(Na,1) * tf - (2 * R_ref / C ./ D)' * ones(1,Nr); 
Ssc = exp(1i * pi * (Km .* alpha)'* ones(1,Nr) .* tau.^2);  %CS变标方程

Sr_rA = ftx(Sr);                       % 变换到距离多普勒域  第一步
Sr_sc = Sr_rA .* Ssc;                  %乘以变标方程        第二步

Sr_RA = fty(Sr_sc);                    %变换到二维频域      第三步
rangeMod_phase = pi * (D ./ Km)' * ones(1,Nr) .* (ones(Na,1) * fr.^2);  %距离调制相位   
cs_bulk_phase = 4 * pi * R_ref / C * alpha' * ones(1,Nr) .* (ones(Na,1) * fr); %一致RCM相位
Sr_RA_cor = Sr_RA .* exp(1i * (rangeMod_phase + cs_bulk_phase));  %相位相乘    第四步

Sr_rA_1 = ifty(Sr_RA_cor);   %变换到距离多普勒域      第五步
dop_phase = 4 * pi / lambda * (ones(Na,1) * R) .* (D' * ones(1,Nr)); %方位压缩
phase_cor = 4 * pi / C^2 * (Km .* ((1 - D)./D .^ 2))' * ones(1,Nr) .* (ones(Na,1) * (R_ref - R) .^ 2); %附加相位校正
Sr_rA_cor = Sr_rA_1 .* exp(1i * (dop_phase + phase_cor));   % 相位相乘  第六步
% Sr_rA_cor = Sr_rA_1 .* exp(1i * dop_phase);   % 相位相乘  第六步
Sr_ra = iftx(Sr_rA_cor);  %变换到时域  第七步

%========================画图============================================

%================显示场景点====================================
figure;
subplot(211);
plot(Targets(:,2),Targets(:,1),'rO');
grid on;
axis([Scene_center(2) + Y_min Scene_center(2) + Y_max Scene_center(1) + X_min Scene_center(1) + X_max]);
xlabel('距离向－AT （m）');
ylabel('方位向－CT （m）');
title('目标场景坐标');

% figure;
% imagesc(abs(Sr_rA));
% title('距离多普勒域-原始信号');
% 
% figure;
% imagesc(abs(Sr_sc));
% title('距离多普勒域-变标后信号');
% 
% figure;
% imagesc(abs(Sr_RA));
% title('二维频域-变标后信号');
% 
% figure;
% imagesc(abs(Sr_RA_cor));
% title('二维频域-补余RCM后');
% 
% figure;
% imagesc(abs(Sr_rA_1));
% title('距离多普勒域-补余RCM后');
% 
% figure;
% imagesc(abs(Sr_rA_cor));
% title('距离多普勒域-方位压缩和相位校正后');
figure;
mesh(abs(Sr_ra));

figure;
imagesc(abs(Sr_ra));
title('时域-点目标');
