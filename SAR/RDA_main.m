% 2016/10/17
%Liu Yakun 

clc;
clear;
close all;

% 雷达平台参数
C = 3e8; %光速
Fc = 5e9;%载频
lambda = C / Fc;%波长
Vr = 150;%飞行速度
%theta_rc = 0 / 180 * pi; %波束斜视角
H = 5000;%飞行高度
La = 4;
%R0 = 2e4;%


%场景参数
Y0 = 10000;%场景中心Y轴
R0 = sqrt(H^2 + Y0^2);%航线到场景中心最近斜距
Length_X = 300;
Length_Y = 800;
Beta = atan(Y0 / H);% 场景中心入射角
scene_center = [0,Y0];

theta = lambda / La;%波束宽度
Lsar = theta * R0;%合成孔径长度
Tsar = Lsar / Vr;%合成孔径时间

% 快时间参数
Tr = 2.5e-6;
Kr = 20e12;
alpha_Fsr = 1.2;%  距离过采样率
Fs_org = alpha_Fsr * Kr * Tr;%距离向原始采样率
Ts_org = 1 / Fs_org;
Rmin = sqrt((Y0 - Length_Y/2)^2 + H^2);%最近斜距
Rmax = sqrt((Y0 + Length_Y/2)^2 + H^2 + (Lsar/2)^2);%最远斜距
sample_time = 2 * (Rmax - Rmin) / C + Tr;%距离采样时间长度
Nr_org = sample_time / Ts_org;% 距离向采样点数 
Nr = 2^nextpow2(Nr_org);
Ts = sample_time / Nr;% 更新后的采样间隔
% Tf_org = [-Nr / 2:(Nr / 2 -1)] * Ts;%中心为0的采样时间向量
% Tf = Tf_org + 2 * R0 / C;% 距离向采样时间矩阵
Tf = linspace(2 * Rmin / C,2 * Rmax / C + Tr,Nr);%快时间采样时间向量
Rf = Tf * C / 2;%斜距向量

% 慢时间参数
alpha_Fsa = 1.25;
Ka = -2 * Vr^2 / lambda /R0;%方位向条频率
PRF_org = alpha_Fsa * Ka * Tsar;%原始PRF
Na_org = (Length_X + Lsar) / Vr * PRF_org;%方位向采样数
Na = 2^nextpow2(Na_org);%为了做FFT 更新的
PRF = (Length_X + Lsar) / Vr / Na;%最终的PRF
Ts = linspace(-Lsar / 2 / Vr,(Length_X + Lsar / 2) / Vr,Na);
Ra = Ts * Vr;

Targets = [ 0  -100   1
            200 0     1
            0   100   1
           -200  0    1];
nTargets = size(Targets,1);
Targets(:,1:2) = Targets(:,1:2) + ones(nTargets,1) * scene_center;

%脉冲接收
echo = echo_creation(C,H,Y0,lambda,Lsar,Kr,Tr,Tf,Ra,Targets);

x = Y0 + (Tf * C / 2 - R0) / sin(Beta);
y = Ra;
% mesh(abs(echo));

%距离压缩
t = Tf - 2 * Rmin / C; 
signal_ref = exp(1i * Kr * t.^2) .* (t > 0 & t < Tr);
signal_compressed = pulse_compression(echo,signal_ref,Na);

% mesh(abs(signal_compressed));
colormap(gray);
imagesc(x,y,255-abs(signal_compressed));
