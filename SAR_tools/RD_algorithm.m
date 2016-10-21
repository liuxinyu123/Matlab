%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      SAR imaging: RD_algorithm (2010.08)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
close all;
%%%%%% 参数设置
%%%%%%%%%%%%%%% constant parameters %%%%%%%%%%%%%%%%%%%

c = 3e8;                   % speed of light,3*10^8m/s

%%%%%%%%%%%%%%% chirp signal parameters %%%%%%%%%%%%%%%

fc = 1e9;                % carrier frequency,1GHz
wc = 2 * pi * fc;
lambda = c / fc;           % Wavelength at carrier frequency
Tr = 1.5e-6;               % chirp pusle duration 1.5us
Br = 150e6;                % chirp frequebcy bandwidth 150Mhz
Kr = Br/Tr;                % chirp signal: frequency modulation rate

%%%%%%%%%%%%%% observating strip parameters %%%%%%%%%%%%

H = 3000;                          % 飞行平台的高度
ThetaCenter = 45 / 180 * pi;       % 正侧视SAR系统的斜视角？？？ 擦地角的余角
Platform_center = [0, 0, H];       % 飞行平台的位置坐标

Y_C = H * tan(ThetaCenter);
Scene_center = [0, Y_C, 0];        % 场景中心的位置坐标
delta_X = 150;                     % 场景区域的范围：x坐标（平台飞行方向）
delta_Y = 100;                     % 场景区域的范围：y坐标
Xmin = -delta_X / 2;
Xmax = delta_X / 2;
Ymin = -delta_Y / 2;
Ymax = delta_Y / 2;
RC0 = sqrt(sum((Scene_center - Platform_center) .^2));  %天线平台中心到场景中心的距离

%%%%%%%%%%%%%% performance parameters %%%%%%%%%%%%%%%%%%

rho_R = c / (2 * Br);             % 距离分辨率（m）
roh_Y = rho_R / sin(ThetaCenter); % 透视矫正后的切航迹向(CT)分辨率（m）
rho_AT = rho_R;                   % 沿航迹向(AT)的理论分辨率（m）

V = 150;                         % 飞行平台的速度（x轴：m/s）
D = 2 * rho_AT;                  % 沿航迹向的实天线孔径长度，page 63, function(3.27)
Lsar = lambda * RC0 / D;         % AT向一个合成孔径的长度; page 59, function(3.14)
Tsar = Lsar / V;                 % 一个合成孔径时间

%%%%% fast-time sampling sequence

Rate_Fs = 1.2;                   % 快时间域的过采样系数为1.2
Fs = round(Rate_Fs * Br);        % 距离维的过采样率（双通道：实部＝I， 虚部＝Q）
Ts = 1 / Fs;                     % 快时间采样时间间隔
delta_Rs = Ts * c;               % 快时间采样间隔对应的距离长度

Rmin = sqrt((Y_C + Ymin)^2 + H^2);                  % 场景最短距离
Rmax = sqrt((Y_C + Ymax)^2 + (Lsar/2)^2 + H^2);     % 场景最大距离
Nfast = ceil((2 * (Rmax-Rmin) / c + Tr) / Ts);      % 距离向的采样点个数
Nf = 2^nextpow2(Nfast);                             % for fft
tf_org = [-Nf/2 : (Nf/2 - 1)] * Ts;                 % 理想快时间采样序列
% tff = tf_org*c;
% tr_ff = (tff/2 + RC0) / sin(pi/4); 
tf = (2 * RC0 / c) + tf_org;                        % 实际快时间采样值
tr = tf * c / 2;                                    % 快时间采样对应的距离域(单程距)

%%%%% slow-time sampling sequence

Ka = -2 * V^2 / (lambda * RC0);  % 多普勒调频率,(美SAR，p159, func(6.5))
Ba = abs(Ka * Tsar);             % 多普勒带宽
Rate_PRF = 1.25;                 % 慢时间域的过采样系数为1.25
PRF = round(Rate_PRF * Ba);      % 脉冲重复频率
PRT = 1 / PRF;                   % 脉冲重复周期,（一个脉冲重复周期内，飞机的飞行时间间隔）

Nslow = ceil((delta_X + Lsar) / V / PRT);            % 方位向采样点个数
Ns  = 2^nextpow2(Nslow);                             % 优化采样点数（for fft）
% Ns  = Nslow;
ts = (-Ns/2 : (Ns/2 - 1)) * PRT;                     % 慢时间采样序列
ta = ts * V;                                         % 慢时间采样对应的距离域

%%%%%
% % delta_d_AT = lambda;           % AT向的阵元间隔（由SAR原理形成）
% % PRF = round(V / delta_d_AT);   % 脉冲重复频率
% % p = nextpow2(A) returns the smallest power of two that is greater than or equal to the absolute value of A.
% X_start = Xmin - Lsar/2;                             % 场景范围：x轴最小值
% X_end = Xmax + Lsar/2;                               % 场景范围：x轴最大值
% ts = linspace(X_start / V, X_end / V, Ns);           % 对方位向时间采样
% % temp_ts = diff(ts);
% % PRT = sum(temp_ts(:)) / length(temp_ts);               % 离散化脉冲重复周期
% dts = ts(2) - ts(1);              % 离散化脉冲重复周期
% Fsa = 1/dts;                      % 脉冲重复周期
% delta_d_AT = V * Fsa;             % AT向的阵元间隔

%%%%%%%%%%%%% SAR　Resolution  %%%%%%%%%%%%%%%%%%%%%%%%%

% Dr=c/2/Br;                      % range resolution
% Dx=v/Ba;                        % cross-range resolution

%%%%%%%%%%%% set point targets parameters %%%%%%%%%%%%%%%

% format [x, y, z, RCS]          
Ptargets = [  0,   0,   0,   1/2  %];               
            -20,   0    0,   1/2
             50,   0    0,   1/2
              0,  20,   0,   1/2
             20,  20,   0,   1/2];   % 定义点目标的相对位置及兵役散射系数

Ntar = size(Ptargets, 1);          % 点目标个数
Ptargets(:, 1:3) = Ptargets(:, 1 : 3) + ones(Ntar, 1) * Scene_center;  % 点目标的真实坐标位置

%%%%%%%%%%%%%%%%%%%  Show the targets  %%%%%%%%%%%%%%%%%%%

figure;
plot(Ptargets(:, 1), Ptargets(:, 2), 'o','MarkerEdgeColor','b','MarkerFaceColor','g'); %, 'MarkerSize', 12
grid on
% axis equal;
axis([Scene_center(1)+Xmin Scene_center(1)+Xmax Scene_center(2)+Ymin Scene_center(2)+Ymax]);
xlabel('x轴(AT向)'); 
ylabel('y轴(CT向)');
title('观测场景及点目标');

%=========================================================================%
%=========================================================================%

%%%%%%%%%%%%%%%%%%  Show the launch signal %%%%%%%%%%%%%%%%%
% figure;
% duty = Tr/PRT*100;
% A0 = 1;
% t = 0:Ts:4*PRT;
% rect = A0/2*square(2*pi*PRF*t,duty)+A0/2;
% subplot(2,2,1);
% plot(t,rect);
% grid on;
% axis([0 4*PRT -0.1 A0+0.2]);
% xlabel('x轴-时间');
% ylabel('y轴-幅度');
% title('矩形波（4个周期）');
% 
% t = 0:Ts:2*Tr;
% signal_car = cos(2*pi*fc*t);
% subplot(222);
% plot(t,signal_car);
% grid on;
% axis([0 2*Tr -1.2 1.2]);
% xlabel('x轴-时间');
% ylabel('y轴-幅度');
% title('载波（局部 时间2个脉冲持续时间）');
% 
% t = 0:Ts:2*Tr;
% signal_mod = cos(2*pi*fc*t+Kr*t.*t);
% subplot(223);
% plot(t,signal_mod);
% grid on;
% axis([0 2*Tr -1.2 1.2]);
% xlabel('x轴-时间');
% ylabel('y轴-幅度');
% title('载波+线性调频波（局部 时间为2个脉冲持续时间）');
% 
% signal_launch = rect(1:size(signal_mod,2)).*signal_mod;
% subplot(224);
% plot(t,signal_launch);
% grid on;
% axis([0 2*Tr -1.2 1.2]);
% xlabel('x轴-时间');
% ylabel('y轴-幅度');
% title('发射信号（局部 时间为2个脉冲持续时间）');


%=========================================================================%
%=========================================================================%
disp('*********************************************************');
fprintf(1,'系统参数显示界面：\n');
fprintf(1,'\n');
disp(['快时间域采样点数是:', poly2str(Nf)]);
disp(['快时间域的过采样率是:', poly2str(Rate_Fs)]);

disp(['PRF:', poly2str(PRF)]);
disp(['慢时间域采样点数是:', poly2str(Ns)]);
disp(['慢时间域的过采样率是:', poly2str(Rate_PRF)]);

disp(['AT向合成孔径的长度是:', poly2str(Lsar), ' m']);
disp(['距离向的分辨率是:', poly2str(rho_R), ' m']);
disp(['沿航向的分辨率是:', poly2str(rho_AT), ' m']);
fprintf(1,'\n');
disp('*********************************************************');
          
%=========================================================================%
%=========================================================================%

%%%%%%%%%%%%%%  Generate the raw signal data  %%%%%%%%%%%%%%%%%%%

Echo_data = zeros(Ns, Nf);                  % 初始化回波数据域
for ii = 1 : Ntar
    sigma = Ptargets(ii, 4);                % 获取反射系数
    Xslow = ts .* V - Ptargets(ii, 1);                                % 雷达在慢时域的移动距离，Xslow为1×Ns的矩阵
    Yslow = Platform_center(2) - Ptargets(ii, 2);
    Zslow = Platform_center(3) - Ptargets(ii, 3);
    R = sqrt(Xslow.^2 + Yslow^2 + Zslow^2); % 雷达与目标的距离
    tau = 2 * R / c;                        % 信号走双程的总延时，tau为1*N的矩阵
    Dfast = ones(Ns,1) * tf - tau' * ones(1, Nf);                      % 雷达相对目标移动在快时域产生的时间差，tm为1*M的矩阵，Dfast就为N*M的矩阵
    phase = pi*Kr*Dfast.^2 - (4 * pi / lambda * R') * ones(1,Nf);     % 第一项是指信号走双程的总延时所带来的相位延迟，第二项指双程距离R所带来的相位延迟
    Echo_data = Echo_data + sigma * exp(1i*phase) .* (abs(Dfast) <= Tr/2) .* ((abs(Xslow) <= Lsar/2)' * ones(1,Nf));
    
    %%%%%%%%%%%%%%%%%%%%% Show the perproity about targets %%%%%%%%%%%%%%%%%%%%%%%
    
%     figure;
%     subplot(221);
%     plot(ts,abs(Xslow));
%     grid on;
%     xlabel('X轴－时间');
%     ylabel('Y轴－值');
%     title(['X轴距离差-目标 ' num2str(ii)]);
%     
%     subplot(222);
%     plot(ts,R);
%     grid on;
%     xlabel('X轴－时间');
%     ylabel('Y轴－值');
%     title(['雷达与目标距离-目标 ' num2str(ii)]);
%     
%     subplot(223);
%     plot(tr,tf);
%     grid on;
%     xlabel('X轴－时间');
%     ylabel('Y轴－值');
%     title(['快时间采样时间-目标 ' num2str(ii)]);
%     
%     subplot(224);
%     plot(ts,tau);
%     grid on;
%     xlabel('X轴－时间');
%     ylabel('Y轴－值');
%     title(['慢时间时间差-目标 ' num2str(ii)]);
%     
%     figure;
%     xx = ta;
%     yy = Y_C + (tr-RC0) ./ sin(ThetaCenter);
%     [x,y] = meshgrid(xx,yy);
%     subplot(221);
%     mesh(x,y,tau' * ones(1, Nf));
%     xlabel('X轴－方位向');
%     ylabel('Y轴－距离向');
%     title(['慢时间延迟-目标 ' num2str(ii)]);
%     
%     subplot(222);
%     mesh(x,y,ones(Ns,1) * tf);
%     xlabel('X轴－方位向');
%     ylabel('Y轴－距离向');
%     title(['时间矩阵-目标 ' num2str(ii)]);
%     
%     subplot(223);
%     mesh(x,y,phase);
%     xlabel('X轴－方位向');
%     ylabel('Y轴－距离向');
%     title(['相位-目标 ' num2str(ii)]);
%     
%     subplot(224);
%     mesh(x,y,abs(Echo_data));
%     xlabel('X轴－方位向');
%     ylabel('Y轴－距离向');
%     title(['回波-目标 ' num2str(ii)]);
    
  
    
    
    
end

%%%%%%%%%%%%%%%%%   Range compression   %%%%%%%%%%%%%%%%%%%%%%

%%%%% 窗函数
% WinTr = (abs(tf_org) <= Tr/2);
% WinIndex = find(WinTr ~= 0);
% WinLength = length(WinIndex);
% WinStart = min(WinIndex(:));
% WinEnd = max(WinIndex(:));
% WindowTr = zeros(size(tf));
% WindowTr(1, WinStart : WinEnd) = hamming(WinLength).';
% % figure,plot(WindowTr);
% h_ref = exp(j * pi * Kr * tf_org.^2) .* WindowTr;                     % 距离向参考函数,时域

h_ref = exp(1i * pi * Kr * tf_org.^2) .* (abs(tf_org) <= Tr/2);          % 距离向参考函数,时域
h_ref = h_ref .* (hamming(Nf).');                                       % 加窗后的参考函数
H_ref = fty(ones(Ns, 1) * h_ref);
% f = linspace(-Fs/2,Fs/2,Nf);
% H_ref = ones(Ns, 1) * exp(1i*pi*f.^2/Kr);
figure, plot(abs(H_ref(2, :)));
title('参考函数的频谱');

Comp_f = fty(Echo_data) .* conj(H_ref);                                 % 距离压缩的频域形式
Comp_tsf = ifty(Comp_f);                                                % 距离压缩的“快时间-慢时间”域形式（ts-tf）
Comp_Rfd = ftx(Comp_tsf);                                               % Azimuth FFT and Range-Doppler domain (fd-tf)
figure, plot(abs(Comp_tsf(Ns/2, :)));
title('距离压缩后的数据抽检');

%%%%%%%%%%%%%   距离弯曲（徙动）校正 -- RCMC   %%%%%%%%%%%%%%%%%%

fd_r = [-Nf/2 : (Nf/2 - 1)] * Fs / Nf;
FF = ones(Ns, 1) * fd_r;                                 % FF为N*M的矩阵
fdc = 0;                                                 % doppler center
fd_a = [-Ns/2 : (Ns/2 - 1)] * PRF / Ns;
FU = fd_a.' * ones(1, Nf);
Refcorr = exp(1i * pi / fc^2 / Ka * (FU.*FF).^2 + 1i * pi * fdc^2 / fc / Ka * FF - 1i * pi / fc / Ka * FU.^2 .* FF); % Range-Doppler domain

%RCMC function

RanComff=ftx(Comp_f);                 % FFT in Azimuth 
RanComffcorr=RanComff .* Refcorr;     % RCMC
RanComtfcorr=ifty(RanComffcorr);      % data in Range_Doppler domain
RanComttcorr=iftx(RanComtfcorr);      % data in Range_Azimuth domain

%%Azimuth compression
ts_mid = ts - 0/V;                                                          % 与场景参考点的慢时域时间差
Refa = exp(1i * pi * Ka * ts_mid.^2) .* (abs(ts_mid) < Tsar/2);              % 方位压缩参考函数
Sa = iftx(Comp_Rfd  .* (conj(ftx(Refa)).' * ones(1, Nf)));                  % 未进行距离徙动校正后的距离方位压缩结果
Sa_RCMC = iftx(ftx(RanComttcorr).*(conj(ftx(Refa)).'*ones(1, Nf)));         % 对距离徙动校正后的距离方位压缩结果



%%%%%%%%%%%%%%%%%%%   成像过程中的图像显示   %%%%%%%%%%%%%%%%%%%%%
% %%%% 01 原始回波数据
figure,
G = 20 * log10(abs(Echo_data));              % 换算成分贝（dB）的形式显示
gm = max(max(G));
thr01 = 40;                                  % 显示的动态范围40dB,"thr" = threshold
g_thr = gm - thr01;
G = (G - g_thr) * (255 / thr01) .* (G > g_thr);
imagesc(tr, ta, -G);                         % 显示原始回波数据图像 
colormap(gray);                              % 使显示的图像为灰度图像             
grid on,axis tight,
xlabel('Range')
ylabel('Azimuth')
title(['(a)原始信号, 场景中心的单程距为：Rc = ', num2str(RC0), 'm'])
% 
%%%% 02 距离压缩后数据（未进行 RCMC 前）
figure;
G = 20 * log10(abs(Comp_tsf));              % 换算成分贝（dB）的形式显示
gm = max(max(G));
thr02 = 40;                                  % 显示的动态范围40dB,"thr" = threshold
g_thr = gm - thr02;
G = (G - g_thr) * (255 / thr02) .* (G > g_thr);
imagesc(tr, ta, -G);
colormap(gray);                              % 使显示的图像为灰度图像             
grid on,axis tight,
xlabel('Range')
ylabel('Azimuth')
title(['(b)距离压缩后的时域信号, Rc = ', num2str(RC0), 'm'])

%%%% 03 距离压缩后距离多普勒域（未进行 RCMC 前）
figure,
G = 20 * log10(abs(Comp_Rfd));               % 换算成分贝（dB）的形式显示
gm = max(max(G));
thr03 = 40;                                  % 显示的动态范围40dB,"thr" = threshold
g_thr = gm - thr03;
G = (G - g_thr) * (255 / thr03) .* (G > g_thr);
imagesc(tr, fd_a, -G);
colormap(gray);                              % 使显示的图像为灰度图像             
grid on,axis tight,
xlabel('Range')
ylabel('Doppler')
title('RCMC前：Range Doppler domain')

%%%% 04 RCMC后: 距离压缩后距离多普勒域
figure,
G = 20 * log10(abs(RanComtfcorr));           % 换算成分贝（dB）的形式显示
gm = max(max(G));
thr04 = 40;                                  % 显示的动态范围40dB,"thr" = threshold
g_thr = gm - thr04;
G = (G - g_thr) * (255 / thr04) .* (G > g_thr);
imagesc(tr, fd_a, -G); 
colormap(gray);                              % 使显示的图像为灰度图像             
grid on,axis tight,
xlabel('Range')
ylabel('Doppler')
title('RCMC后：Range Doppler domain')

%%%% 05 未进行距离徙动校正后的距离方位压缩结果
figure,
G = 20 * log10(abs(Sa));                     % 换算成分贝（dB）的形式显示
gm = max(max(G));
thr05 = 40;                                  % 显示的动态范围40dB,"thr" = threshold
g_thr = gm - thr05;
G = (G - g_thr) * (255 / thr05) .* (G > g_thr);
imagesc(tr, ta, -G);
colormap(gray);                              % 使显示的图像为灰度图像             
grid on,axis tight,
xlabel('Range')
ylabel('Doppler')
title('05 未进行距离徙动校正后的距离方位压缩结果')

%%%% 06 未进行距离徙动校正后的距离方位压缩结果
figure,
G = 20 * log10(abs(Sa_RCMC));                % 换算成分贝（dB）的形式显示
gm = max(max(G));
thr06 = 40;                                  % 显示的动态范围40dB,"thr" = threshold
g_thr = gm - thr06;
G = (G - g_thr) * (255 / thr06) .* (G > g_thr);
imagesc(tr, ta, -G); 
colormap(gray);                              % 使显示的图像为灰度图像             
grid on,axis tight,
xlabel('Range')
ylabel('Doppler')
title(['06 对距离徙动校正后的距离方位压缩结果, Rc = ', num2str(RC0), 'm']);

%%%% 07 透视校正后的最终结果
figure,
% tr_rectify = sqrt((tr).^2 - H.^2);
% tr_rectify = Y_C + (tr-RC0) ./ sin(acos(H./(tr)));
tr_rectify = Y_C + (tr-RC0) ./ sin(ThetaCenter);
imagesc(ta, tr_rectify, -G.'); 
% imagesc( -G.'); 
colormap(gray);                              % 使显示的图像为灰度图像             
grid on,axis tight,
xlabel('x轴：AT向')
ylabel('y轴：CT向')
title(['07 透视校正后的最终结果, Y0 = ', num2str(Y_C), 'm']);

