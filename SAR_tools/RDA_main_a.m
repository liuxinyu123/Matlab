% 2016/10/20
%Liu Yakun
clc;clear;close all;
%雷达参数,平台参数与调频信号参数
C = 3e8; %光速
Fc = 1e9;%发射信号载频
lambda = C / Fc;%发射信号波长
La = 4;%雷达天线方位向长度
H = 5000;%雷达平台飞行高度
V = 150;%雷达飞行速度
Theta_center = 45 / 180 * pi;%波束中心入射角
%----------------------------

R0 = H / sin(Theta_center);%波束中心最短斜距
wave_angle = lambda / La;%方位向波束宽度（弧度）
Lsar = wave_angle * R0;%合成孔径长度
Ka = -2*V^2/lambda/R0;%方位向调频率
Tsar = Lsar/V;%合成孔径时间

%场景参数
Xmin = -50;%
Xmax = 50;
Ymin = -300;
Ymax = 300;

%--------------------------------------------
Yc = sin(Theta_center) * R0;%场景中心y轴坐标
Rmin = sqrt(H^2 + (Yc + Ymin)^2);%场景最短斜距
Rmax = sqrt(H^2 + (Yc + Ymax)^2 + (Lsar/2)^2);%场景最大斜距
scene_center = [0,Yc];%
%快时间参数
Tr = 2.5e-6;%脉冲持续时间
Kr = 2e13;%距离向调频率
Fr_alpha = 1.2;%距离向过采样系数

%-----------------------------------------------
Br = Kr * Tr;%脉冲带宽
Fr = Fr_alpha * Br;%距离向采样频率
Nr = ceil(((Rmax - Rmin) * 2 / C + Tr) * Fr);%距离向采样点数
Nr = 2^nextpow2(Nr);%为了做fft 
Tfast = linspace(2*Rmin/C,2*Rmax/C+Tr,Nr);%距离向采样时间向量
TFr = (2*(Rmax-Rmin)/C+Tr)/Nr;
% Fr = 1/TFr;
Rr = Tfast * C / 2;%距离向斜距矩阵


%慢时间参数
Fa_alpha = 1.3;%方位向过采样系数

%------------------------------------------------
Fa = Fa_alpha * abs(Ka * Tsar);%方位向采样频率
PRF = Fa;%脉冲重复发射周期
Na = ceil(((Xmax - Xmin) + Lsar)/V*PRF);%方位向采样点数
Na = 2^nextpow2(Na);
Tslow = linspace((Xmin-Lsar/2)/V,(Xmax+Lsar/2)/V,Na);%方位向采样时间矩阵
Ra = Tslow * V;%方位向采样对应的距离矩阵
TFa = ((Xmax-Xmin)+Lsar)/V/Na;
% PRF = 1 / TFa;


%点目标坐标
Targets = [ 20   -50    1
            20    50    1
             0    0     1
           -20    50    1
           -20   -50    1]; %格式为场景中（地距） 距离向坐标 方位向坐标 RCS
      
%------------------------------------------
Targets(:,1:2) = Targets(:,1:2) + ones(size(Targets,1),1) * scene_center;
DX = La/2;%方位向分辨率
DR = C/2/Br;%距离向分辨率
DY = DR / sin(Theta_center);

%回波产生
N = size(Targets,1);%目标数量
echo = zeros(Na,Nr); %初始化回波数据
for i = 1:N
    delta_x = Ra - Targets(i,1); %平台距目标的x轴距离差
    delta_y = Targets(i,2);%平台距目标的y轴距离差
    range = sqrt(delta_x.^2 + H^2 + delta_y^2); % 每个方位采样点 雷达到目标的斜距
    rcs = Targets(i,3); %点目标的后向散射系数
    tau = ones(Na,1) * Tfast - (2*range / C)' * ones(1,Nr); % 时间差矩阵
    phase = pi * Kr * tau.^2 - 4 * pi / lambda * (range' * ones(1,Nr));%接收信号相位
    echo = echo + rcs * exp(1i * phase) .* ((abs(delta_x) < Lsar / 2)' * ones(1,Nr)) .* (tau > 0 & tau < Tr);
end

%距离压缩
t = Tfast - 2*Rmin/C;
refr = exp(1i*pi*Kr*t.^2) .* (t > 0 & t < Tr);
signal_comp = ifty(fty(echo) .* (ones(Na,1)*conj(fty(refr))));%参考信号补零后DFT 取共轭



%距离徙动校正
signal_rD = ftx(signal_comp);
win = waitbar(0,'最近邻域插值');
for i = 1:Na
    for j = 1:Nr
        delta_R = (1/8)*(lambda/V)^2*(R0+(j-Nr/2)*C/2/Fr)*((j-Nr/2)/Nr*PRF)^2;%距离徙动量
        RCM = delta_R / DY;%徙动了多少个距离单元
        delta_RCM = RCM - floor(RCM);%小数部分
        if round(RCM + j) > Nr
            signal_rD(i,j) = signal_rD(i,Nr/2);
        else
            if delta_RCM < 0.5
                signal_rD(i,j) = signal_rD(i,j+floor(RCM));
            else
                signal_rD(i,j) = signal_rD(i,j+ceil(RCM));
            end
        end
    end
        waitbar(i/Na);
end
close(win);
signal_rcm = iftx(signal_rD);   


%方位向压缩
ta = Tslow;
refa = exp(1i*pi*Ka*ta.^2).*(abs(ta) < Tsar/2);
final_signal = iftx(ftx(signal_rcm) .* (conj(ftx(refa)).' * ones(1,Nr)));
%绘图

%点目标坐标
figure;
plot(Targets(:,2),Targets(:,1),'ro');
grid on;
axis([Yc+Ymin Yc+Ymax Xmin Xmax]);
xlabel('距离向（米）');
ylabel('方位向（米）');
title('场景中点目标位置');

figure;
subplot(211);
colormap(gray);
imagesc(Rr,Ra,255-abs(signal_comp));
xlabel('距离向（米）');
ylabel('方位向（米）');
title('脉冲压缩后的信号');
subplot(212);
colormap(gray);
imagesc(Rr,Ra,255-abs(signal_rcm));
xlabel('距离向（米）');
ylabel('方位向（米）');
title('RCMC后的信号');



figure;
colormap(gray);
xx = Yc + (Rr-R0)/sin(Theta_center);
yy = Ra;
imagesc(xx,yy,255-abs(final_signal));
xlabel('距离向（米）');
ylabel('方位向（米）');
title('最终的点目标');

