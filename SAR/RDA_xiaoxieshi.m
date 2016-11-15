%2016／10／31
%Liu Yakun
clc,clear,close all;

C = 3e8;%光速
Fc = 1e9; %载波频率
lambda = C / Fc;%波长
H = 5000;%雷达高度
Theta = 2 / 180 * pi;%斜视角
V = 150;%雷达速度
La = 4;%雷达方位向天线长度
wave_width = lambda / La;%波束宽度
Tr = 5e-6;%脉冲持续时间
Br = 60e6;%脉冲带宽
Kr = Br / Tr;%距离向调频率

Y0 = 5e3;%场景中心Y轴距离
Ymin = Y0 -100;
Ymax = Y0 + 100;
Xmin = 0;
Xmax = 50;

R0 = sqrt(H^2 + Y0^2);%最短斜距
Rc = R0 / cos(Theta);%波束中心穿越时刻斜距
Lsar = wave_width * Rc / cos(Theta);%合成孔径长度
Tsar = Lsar / V;%合成孔径时间

Targets = [20   -20    1
           20    20    1
           40    -40   1
           40    40    1];
scene_center = [0   Y0];%场景参考位置
Targets(:,1:2) = Targets(:,1:2) + ones(size(Targets,1),1) * scene_center;

Rmin = sqrt(H^2 + Ymin^2);%场景距雷达最近距离
Rmax = sqrt(H^2 + Ymax^2 + (Lsar/2)^2);%场景距雷达最远距离
% alpha_fsr = 1.2;%距离向过采样系数
% alpha_fsa = 1.3;%方位向过采样系数
fs = Br;%距离向采样频率
Nr = (2 * (Rmax - Rmin) / C + Tr) * fs;%距离向采样点数
Nr = 2^nextpow2(Nr);
fs = Nr / (2 * (Rmax - Rmin) / C + Tr);
alpha_fsr = fs / Br;
% alpha_fsr = Nr / (2 * (Rmax - Rmin) / C + Tr) / Br;
tf = linspace(2 * Rmin / C,2 * Rmax / C + Tr,Nr);%距离向快时间向量
Ts = 1 / fs;

Ka = -2 * V^2 * cos(Theta)^3 / lambda / R0;%方位向调频率
Ba = round(abs(Ka * Tsar));%方位向多普勒频率
PRF = Ba;%脉冲重复频率
Na = (Xmax - Xmin + Lsar) / V * PRF;%方位向采样点数
Na = 2^nextpow2(Na);
alpha_fsa = Na / 
ts = linspace((Xmin - R0 * tan(Theta) - Lsar / 2)/V,(Xmax + Lsar/2)/V,Na);%方位向慢时间向量

Ntar = size(Targets,1);
echo = zeros(Na,Nr);

for i = 1:Ntar
    rcs = Targets(i,3);
    delta_x = abs(ts*V - Targets(i,1));
    delta_y = Targets(i,2);
    delta_z = H;
    
    R = sqrt(delta_x.^2 + delta_y^2 + delta_z^2);%瞬时斜距
    tau = 2 * R / C;%斜距对应的传播时间
    delta_t = ones(Na,1) * tf - tau.' * ones(1,Nr);
    phase = -4*pi/lambda*R.'*ones(1,Nr) + Kr*pi*delta_t.^2;
    echo = echo + rcs * exp(1i*phase) .* (delta_t > 0 & delta_t < Tr) .* ((abs(delta_x)<(Lsar/2 + R0*tan(Theta))).'*ones(1,Nr));
end

% figure;
% colormap(gray);
% imagesc(255-abs(echo));

%距离脉压
t = tf - 2 * Rmin / C;
ref_r = exp(1i*pi*Kr*t.^2) .* (t < Tr & t > 0);
signal_comp = ifty(fty(echo) .* (ones(Na,1) * conj(fty(ref_r))));
figure;
colormap(gray);
imagesc(255-abs(signal_comp));

%RCMC
signal_RD = ftx(signal_comp);
win = waitbar(0,'最近邻域插值');
for i = 1:Na
    for j = 1:Nr
        rcm = 1/8*(lambda/V)^2*(tf(j)*C/2)*(-Ka*ts(i))^2;
        nrcm = 2*rcm/C/Ts;
        delta_nrcm = nrcm - floor(nrcm);
        
        if round(nrcm)+j > Nr
            signal_RD(i,j) = signal_RD(i,Nr/2);
        else
            if delta_nrcm < 0.5
                signal_RD(i,j) = signal_RD(i,j+floor(nrcm));
            else
                signal_RD(i,j) = signal_RD(i,j+ceil(nrcm));
            end
        end
    end
    waitbar(i/Na);
end
close(win);

signal_rcmc = iftx(signal_RD);

figure;
colormap(gray);
imagesc(255-abs(signal_rcmc));