% 2016/10/17
%Liu Yakun 

clc;
clear;
close all;

% 雷达平台参数
C = 3e8; %光速
Fc = 1e9;%载频
lambda = C / Fc;%波长
Vr = 150;%飞行速度
H = 5000;%飞行高度
La = 4;


%场景参数
Y0 = 10000;%场景中心Y轴
R0 = sqrt(H^2 + Y0^2);%航线到场景中心最近斜距
Length_X = 400;
Length_Y = 800;
Beta = atan(Y0 / H);% 场景中心入射角
scene_center = [Length_X/2,Y0];

theta = lambda / La;%波束宽度
Lsar = theta * R0;%合成孔径长度
Tsar = Lsar / Vr;%合成孔径时间

% 快时间参数
Tr = 2.5e-6;
Kr = 20e12;
alpha_Fsr = 1.2;%  距离过采样率
Fs_org = alpha_Fsr * Kr * Tr;%距离向原始采样率
Tsr_org = 1 / Fs_org;
Rmin = sqrt((Y0 - Length_Y/2)^2 + H^2);%最近斜距
Rmax = sqrt((Y0 + Length_Y/2)^2 + H^2 + (Lsar/2)^2);%最远斜距
sample_time = 2 * (Rmax - Rmin) / C + Tr;%距离采样时间长度
Nr_org = ceil(sample_time / Tsr_org);% 距离向采样点数 
Nr = 2^nextpow2(Nr_org);
% Tsr = sample_time / Nr;% 更新后的采样间隔
% Tf_org = (-Nr / 2:(Nr / 2 -1)) * Tsr_org;%中心为0的采样时间向量
% Tf = Tf_org + 2 * R0 / C;% 距离向采样时间矩阵
Tf = linspace(2 * Rmin / C,2 * Rmax / C + Tr,Nr);%快时间采样时间向量
Rf = Tf * C / 2;%斜距向量

% 慢时间参数
alpha_Fsa = 1.25;
Ka = -2 * Vr^2 / lambda /R0;%方位向条频率
Fdop = abs(Ka * Tsar);%多普勒频率
PRF_org = alpha_Fsa * Fdop;%原始PRF
PRT_org = 1 / PRF_org;
Na_org = ceil((Length_X + Lsar) / Vr / PRT_org);%方位向采样数
Na = 2^nextpow2(Na_org);%为了做FFT 更新的
Tsa = (Length_X + Lsar) / Vr / Na;%方位向采样间隔
% PRF = 1 / Tsa;%最终的PRF
Ts = linspace(-Lsar / 2 / Vr,(Length_X + Lsar / 2) / Vr,Na);
% Ts = (-Na/2:(Na/2-1))*PRT_org;
Ra = Ts * Vr;

Targets = [ 0  -100   1
            100 0     1
            0   100   1
           -100  0    1];
nTargets = size(Targets,1);
Targets(:,1:2) = Targets(:,1:2) + ones(nTargets,1) * scene_center;

DX = La / 2; %方位向分辨率
DY = C / (2 * Kr * Tr); %距离向分辨率
%场景点绘图
figure;
plot(Targets(:,2),Targets(:,1),'ro');
grid on;
xlabel('距离向（米）');
ylabel('方位向(米)');
axis([Y0-Length_Y/2 Y0+Length_Y/2 0 Length_X]);
title('场景点坐标');

%产生回波
echo = echo_creation(C,H,Y0,lambda,Lsar,Kr,Tr,Tf,Ra,Targets);

x = Y0 + (Tf * C / 2 - R0) / sin(Beta);
y = Ra;

%绘制回波数据
% figure;
% mesh(abs(echo));

%距离压缩
t = Tf - 2 * R0 / C; 
refr = exp(1i * Kr * pi * t.^2) .* (abs(t) < Tr/2);
signal_compressed = ifty(fty(echo) .* (ones(Na,1) * conj(fty(refr))));
signal_rD = ftx(signal_compressed);
colormap(gray);
figure;
imagesc(x,y,255-abs(signal_compressed));

%RCMC
signal_RCMC = RCMC(signal_rD,lambda,C,Vr,Tsr_org,DY,PRF_org,2);
figure;
imagesc(x,y,255-abs(signal_RCMC));
title('RCMC后信号');

% 方位向压缩
refa = exp(1i * pi * Ka * Ts.^2) .* (abs(Ts) < Tsar / 2);
final_signal = iftx(ftx(signal_RCMC) .* (conj(ftx(refa)).' * ones(1,Nr)));
figure;
imagesc(x,y,255-abs(final_signal));
title('最终点目标');

% %距离压缩
% tr = Tf - 2*Rmin/C;
% refr = exp(1i*pi*Kr*tr.^2) .* (tr > 0 & tr < Tr);
% signal_comp = ifty(fty(echo) .* (ones(Na,1) * conj(fty(refr))));
% figure;
% imagesc(255-abs(signal_comp));
% 
% %方位压缩
% ta = Ts;
% refa = exp(1i*pi*Ka*ta.^2) .* (abs(ta) < Tsar/2);
% final_signal = iftx(ftx(signal_comp) .* (conj(ftx(refa)).' * ones(1,Nr)));
% figure;
% imagesc(255-abs(final_signal));
