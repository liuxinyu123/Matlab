function [x] = RDA(La,H,Vr,theta_rc,Y0,Fc,Kr,Tr,alpha_Fsr,alpha_Fsa,delta_x,delta_y,Targets)

C = 3e+8;%光速
lambda = C / Fc;%波长

Br = Kr * Tr;%调频信号带宽
Fsr = alpha_Fsr * Br;%距离向采样率
Tsr = 1.0 / Fsr;%采样间隔时间

theta_gc = atan(H / Y0);%场景中心擦地角
R0 = H / sin(theta_gc);%场景中心斜距
Lsar = lambda * R0 / La * cos(theta_rc); %合成孔径长度
Tsar = Lsar / Vr;

Rmin = sqrt(H^2 + (Y0 - delta_y)^2);%最短斜距
Rmax = sqrt(H^2 + (Y0 + delta_y)^2 + (Lsar / 2)^2);%最长斜距
delta_tau = 2 * (Rmax - Rmin) / C;%距离时间间隔
Nr = (delta_tau + Tr) / Tsr;%距离采样点数
Nr_update = 2 ^ nextpow2(Nr);% 为了fft，扩展采样点数到2的幂。
Tf = linspace(2 * Rmin / C,2 * Rmax / C + Tr,Nr_update); %快时间采样时间点
Rf = Tf * C / 2;%快时间采样距离点                   

Ka = 2 * Vr^2 * cos(theta_rc)^3 / lambda / R0; %p93 公式4.38
Ba = Ka * Tsar; %多普勒带宽
Fa = alpha_Fsa * Ba;%方位向采样频率
Tsa = 1.0 / Fa;%方位采样间隔
Na = (2 * delta_x + Lsar) / Vr / Tsa;%方位采样点数
Na_update = 2 ^ nextpow2(Na);%更新后的方位向采样率
Ts = linspace((-delta_x - Lsar / 2) / Vr,(delta_x + Lsar / 2 )/ Vr,Na_update);%方位向采样时间点
Rs = Ts * Vr;%方位向采样对应的距离
R_eta = sqrt(R0^2 + Rs.^2);%瞬时斜距

% t = linspace(-Tr/2,Tr/2,Nr_update);
t = linspace(0,Tr,Nr_update);
replica = exp(1i*pi*Kr*t.^2);%复制信号
figure;
subplot(211);
plot(t,real(replica));
xlabel('时间（秒）');
ylabel('复制信号实部');
grid on;

freq = linspace(-Nr_update/Tr/2,Nr_update/Tr/2,Nr_update);
subplot(212);
plot(freq,abs(fftshift(fft(replica))));
xlabel('频率（Hz）');
ylabel('复制信号频谱');
grid on;

N = size(Targets,1);% 点目标个数
Srev = zeros(Na_update,Nr_update);%初始化接收信号
scene_center = [0 Y0];
Targets(:,1:2) = Targets(:,1:2) + ones(N,1) * scene_center;

figure;
plot(Targets(:,2),Targets(:,1),'ro');
grid on;
xlim([scene_center(2) - delta_y, scene_center(2) + delta_y]);
ylim([-delta_x,delta_x]);
xlabel('距离向（米）');
ylabel('方位向（米)');
title('场景点坐标');

for i = 1:1:N
    sigma = Targets(i,3);%点目标RCS
    Rx = Vr * Ts - Targets(i,1);%点目标与雷达的x轴距离差
    Ry = Targets(i,2);%点目标与雷达的y轴距离差
    R = sqrt(Rx.^2 + Ry^2 + H^2);%点目标的瞬时斜距
    tau = 2 * R / C;%时间延迟
    T_mat = ones(Na_update,1) * Tf - tau' * ones(1,Nr_update);%时间矩阵
    phase = 4*pi*R'/lambda * ones(1,Nr_update) + pi*Kr*T_mat.^2;% 接收信号的相位（正交解调后）
    Srev = Srev + sigma * exp(1i*phase) .* (abs(T_mat) >= 0 & abs(T_mat) <= Tr) .* ((abs(Rx) < Lsar / 2)' * ones(1,Nr_update)) ;
end

row = scene_center(2) + (Rf - R0) / sin(theta_gc);
col = Rs;
figure;
[xx,yy] = meshgrid(row,col);
surf(xx,yy,abs(Srev));

ref_comp = fty(ones(Na_update,1)*(replica .* (t >= 0 & t <= Tr)));%.* (ones(Na_update,1)*hamming(Nr_update)');%复制信号频域加窗
S_comp = fty(Srev) .* conj(ref_comp);%复制信号补零后，FFT之后共轭
S_comp_ra = ifty(S_comp);%  距离压缩后变换到距离时域
S_rD = ftx(S_comp_ra);%变换到距离多普勒域

row = Rf / sin(theta_gc);
col = Rs;
figure;
colormap(gray);
imagesc(row,col,255-abs(S_rD));


