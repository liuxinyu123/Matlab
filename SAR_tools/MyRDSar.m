%=========================================================
% Author: Liu Yakun
% Date: 20,Sep,2016
% Program Name: The image of dot's targets in SAR with RDA
%=========================================================

clc,clear,close all;

%＝＝＝＝＝＝＝＝＝＝＝＝＝参数设置＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
C = 3e+8;       %光速
fc = 1e9;     %发射信号中心频率
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
%             50 -200 0   1
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
    
    %==============接收信号观测=========================================
%     figure;
%     [x,y] = meshgrid(tf * C / 2 / sin(Theta_d),ts * V);
%     surf(x,y,abs(Sr));
%     xlabel('距离向 － CT （m）');
%     ylabel('方位向 － AT （m）');
%     zlabel('接收信号幅度');
%     title('接收信号强度');

end

%====================距离脉压＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
ref_h = exp(1i * pi * Kr * tf_ideal .^ 2) .* (abs(tf_ideal) < Tr / 2);  %脉压参考函数
ref_h = ref_h .* hamming(Nr).';                        %加汉明窗
ref_H = fty(ones(Na,1) * ref_h);                       %变换到频域
Sr_comp_ra = ifty(fty(Sr) .* conj(ref_H));             %频域相乘
Sr_comp_rA = ftx(Sr_comp_ra);                          %方位向傅立叶变换到多普勒域
Sr_RD = Sr_comp_rA;
%查看脉压后的信号
% figure;
% plot(Scene_center(2) + (tr - R0) / sin(Theta_d),abs(Sr_comp_rA(1,:)));
% title('压缩后的数据抽检');

%＝＝＝＝＝＝＝＝＝＝＝＝＝＝距离徙动校正＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
%最近邻插值
% win = waitbar(0,'最近邻域插值');
% 
% for m = 1 : 1 : Na
%     for n = 1 : 1 : Nr
%         delta_R = 1 / 8 * (lambda / V) ^ 2 * (R0 + (n- Nr / 2) * Ts * C / 2) * ((m- Na / 2) / Na * PRF) ^ 2;
%         RMC = 2 * delta_R / C * Fs;
%         delta_RMC = RMC - floor(RMC);
%         
%         if n + round(RMC) > Nr
%             Sr_comp_rA(m,n) = Sr_comp_rA(m,Nr / 2);
%         else
%             if delta_RMC >= 0.5
%                 Sr_comp_rA(m,n) = Sr_comp_rA(m,n + ceil(RMC));
%             else
%                 Sr_comp_rA(m,n) = Sr_comp_rA(m,n + floor(RMC));
%             end
%         end
%     end
%     waitbar(m/Na);
% end
% close(win);

%sinc插值
k = 4;
win = waitbar(0,'sinc插值');
for m = 1 : 1 : Na
    for n = 1 : 1 : Nr
        delta_R = 1 / 8 * (lambda / V) ^ 2 * (R0 + (n- Nr / 2) * Ts * C / 2) * ((m- Na / 2) / Na * PRF) ^ 2;
        RMC = 2 * delta_R / C * Fs;
        delta_RMC = RMC - floor(RMC);
        
        for i = -k / 2 : k / 2 - 1;
            Sr_comp_rA(m,n) = Sr_comp_rA(m,n) + Sr_comp_rA(m,n) .* sinc(n + delta_RMC);
        end;
    end;
    waitbar(m / Na);
end;
close(win);


Sr_RMC = iftx(Sr_comp_rA);

%=====================方位 压缩＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝

ref_a = exp(1i * pi * Ka * ts .^ 2) .* (abs(ts) < Tsar / 2);
Sr_SAR = iftx(ftx(Sr_RMC) .* conj(ftx(ref_a)' * ones(1,Nr)));
% plot(abs(Sr_SAR));

%=======================绘图＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝

%================显示场景点====================================
figure;
subplot(211);
plot(Targets(:,2),Targets(:,1),'rO');
grid on;
axis([Scene_center(2) + Y_min Scene_center(2) + Y_max Scene_center(1) + X_min Scene_center(1) + X_max]);
xlabel('距离向－AT （m）');
ylabel('方位向－CT （m）');
title('目标场景坐标');

%回波信号
colormap(gray);
subplot(212);
row = Scene_center(2) + (tr - R0) / sin(Theta_d);
col = ts * V;
imagesc(row,col,255 - abs(Sr));
xlabel('距离向－AT （m）');
ylabel('方位向－CT （m）');
grid on;
title('回波信号');

%距离脉压之后的数据
figure;
colormap(gray);

subplot(211);
imagesc(row,col,255 - abs(Sr_RD));
xlabel('距离向－AT （m）');
ylabel('方位向－CT （m）');
grid on;
title('距离脉压之后的数据（RCMC之前）');

subplot(212);
imagesc(row,col,255 - abs(Sr_RMC));
xlabel('距离向－AT （m）');
ylabel('方位向－CT （m）');
grid on;
title('距离脉压之后的数据（RCMC之后）');

figure;
colormap(gray);
imagesc(row,col,255 - abs(Sr_SAR));
xlabel('距离向－AT （m）');
ylabel('方位向－CT （m）');
grid on;
title('方位压缩之后的数据');

