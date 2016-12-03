%2016/11/8
%Liu Yakun

clc,clear,close all;
C = 3e8;%光速
% fc = 8e9;%载波频率
% lambda = C / fc;%波长
lambda = 0.03;            %发射信号波长
fc = C / lambda;          %发射信号中心频率
v = 200;%雷达平台速度
% h = 5000;
D = 5;%天线长度
theta = 45 / 180 * pi;%斜视角
beta = lambda / D;%波束宽度
yc = 41.7e3;%场景中心斜距
xc = yc * tan(theta);%以雷达波束中心穿越场景中心点时雷达坐标为（0，0）点，场景中心点的方位向坐标

wa = 100;%方位向宽度
wr = 100;%距离向宽度
xmin = xc - wa/2;%方位向边界点
xmax = xc + wa/2;%方位向边界点
ymin = yc - wr/2;%距离向边界点
ymax = yc + wr/2;%距离向边界点

rnear = ymin/cos(theta-beta/2);%最近斜距
rfar = ymax/cos(theta+beta/2);%最远斜距
rmid = (rnear + rfar) / 2;
a = 0.7;
b = 0.6;

targets = [xc + a*wa/2,yc + b*wr/2
           xc + a*wa/2,yc - b*wr/2
           xc - a*wa/2,yc + b*wr/2
           xc - a*wa/2,yc - b*wr/2
           xc         ,yc         ];


xbegin = xc - wa/2 - ymax*tan(theta+beta/2);%开始照射场景的方位向坐标
xend = xc + wa/2 - ymin*tan(theta - beta/2);%结束照射时的方位向坐标
xmid = (xbegin + xend) / 2;
ka = -2*v^2*cos(theta)^3/lambda/yc;%方位向调频率
% lsar = beta * yc / cos(theta);%合成孔径长
lsar = yc * (tan(theta + beta/2) - tan(theta - beta/2));%合成孔径长
tsar = lsar/v;%合成孔径时间
ba = abs(ka * tsar);%多普勒频率
tr = 2e-6;%脉冲持续时间
br = 50e6;%脉冲带宽
kr = br / tr;%距离向调频率

% PRFmin = ba + 2*v*br*sin(theta)/C;
% PRFmax = 1/(2*tr+2*(rfar-rnear)/C);
% PRF = round(1.3 * ba);
% fs = round(1.2 * br);
alpha_slow = 1.3;%方位向过采样率
alpha_fast = 1.2;%距离向过采样系数
PRF = round(alpha_slow * ba);%重频
PRT = 1 / PRF;%脉冲重复周期
Fs = alpha_fast * br;%距离向采样率
Ts = 1 / Fs;%距离向采样间隔
Na = round((xend - xbegin)/v/PRT);%方位向采样点数
Na = 2^nextpow2(Na);%为了fft，更新点数
Nr = round((tr + 2*(rfar-rnear)/C)/Ts);%距离向采样点数
Nr = 2^nextpow2(Nr);%为了fft，更新点数
% PRF = Na / ((xend - xbegin)/v);
% fs = Nr / (tr + 2*(rfar-rnear)/C);
% tslow = linspace(xbegin/v,xend/v,Na);
ts = [-Na/2:Na/2 - 1]*PRT;%方位采样时间序列
tf = [-Nr/2:Nr/2 - 1]*Ts + 2*rmid/C;%距离采样序列
range = tf*C/2;%距离门




ntargets = size(targets,1);%点目标个数
echo = zeros(Na,Nr);%初始化点目标
for i = 1:ntargets
    xi = targets(i,1);%方位向坐标
    yi = targets(i,2);%距离向坐标
    tci = ((xi - xc) - (yi - yc)*tan(theta)) / v;%波束中心穿越时刻
    rci = yi / cos(theta);%波束中心穿越时刻瞬时斜距
%     tsi = rci*(cos(theta)*tan(theta + beta/2) - sin(theta))/v;%波束开始照射与波束中心照射时刻的方位向时间差
    tsi = yi * (tan(theta + beta/2 - tan(theta))) / v;
%     tei = rci*(sin(theta) - cos(theta)*tan(theta - beta/2))/v;%波束结束照射与波束中心照射时刻的方位向时间差
    tei = yi * (tan(theta) - tan(theta - beta/2)) / v;
    ri = sqrt(yi^2 + (xi - v*ts).^2);%照射时间内的瞬时斜距
    tau = 2 * ri / C;%延时
    t = ones(Na,1)*tf - tau.'*ones(1,Nr);%t-tau矩阵
    phase = pi*kr*t.^2 - 4*pi/lambda*(ri.'*ones(1,Nr));%相位
    
    echo = echo + exp(1i*phase).* (abs(t)<tr/2) .* ((ts > (tci - tsi) & ts < (tci + tei))' * ones(1,Nr));
    
end
 
ff = linspace(-Fs/2,Fs/2,Nr);
% REF_R = exp(-1i * pi * ff.^2 / kr);
% signal_comp = ifty(fty(echo) .* (ones(Na,1) * REF_R));
fdoc = 2*v*sin(theta)/lambda;%多普勒中心频率
% mamb = round(fdoc / PRF);
% fdoc = fdoc - mamb * PRF;
t = tf - 2*rmid/C;%方式二生成匹配滤波器 ：fft之后取共轭
ref_r = exp(1i*pi*kr*t.^2) .* (abs(t) < tr/2);
signal_rfat = fty(echo) .* (ones(Na,1) * conj(fty(ref_r)));%.* (exp(-1i*pi*0.9*fdoc*ts).' * ones(1,Nr));
signal_comp = ifty(signal_rfat);%距离脉压之后的信号
% signal_comp = signal_comp .* (exp(-1i*pi*fdoc*ts).' * ones(1,Nr));
signal_rtaf = ftx(signal_comp);
signal_rfaf = ftx(signal_rfat);%二维信号
% signal_rtat = iftx(ifty(signal_rfaf));
% signal_rtaf = ftx(signal_comp);

% figure;
% subplot(211);
% imagesc(abs(signal_comp));
% subplot(212);
% imagesc(abs(signal_rtat));
fu = linspace(fdoc - PRF/2,fdoc + PRF/2,Na);%方位向频率
d = sqrt(1 - (lambda * fu / 2 / v ).^2);
ksrc = 2*v^2*fc^3*d.^3/C/yc./fu.^2;
H_src = exp(-1i*pi*(ones(Na,1) * ff.^2)./(ksrc.' * ones(1,Nr)));
signal_src = ifty(signal_rfaf .* H_src);
% d = sqrt(1-(lambda*fdoc/2/v)^2);
% ksrc = 2*v^2*fc^3*d^3/C/yc/fdoc;
% ref_R = exp(1i*pi*(1/kr-1/ksrc)*ff.^2);
% signal_comp = ifty(fty(echo) .* (ones(Na,1) * ref_R));

% fu = linspace(fdoc + PRF/2,fdoc - PRF/2,Na);
% fa = fftshift(fu);
% fu = ka * ts;
% d = sqrt(1 - (lambda*fu/2/v).^2);
% ksrc = 2*v^2*fc^3*d.^3/C/yc/fdoc^2;
% km = kr/(1 - kr / ksrc);
% ref_r = exp(1i*pi*km*t.^2) .* (abs(t) < tr/2);
% signal_comp1 = ifty(fty(echo) .* (ones(Na,1) * conj(fty(ref_r))));
% t = tf - 2*rnear/C;
% ref_src = exp(1i*pi*ksrc*t.^2) .* (abs(t) < tr/2);
% H_src = exp(-1i*pi*(ones(Na,1) * ff.^2)./(ksrc.' * ones(1,Nr))); 
% H_src = fty(ref_src);
% H_src = exp(-1i * pi * yc * C / (2 * v^2 * fc^3) * ((fu.^2 ./ d.^3)' * ones(1,Nr)) .* f.^2);
% signal_src = ifty(signal_rfaf .* H_src);

% signal_RD = ftx(signal_comp);%距离多普勒域信号
signal_RD = signal_src;
% signal_rcmc = zeros(Na,Nr);%初始化
% win = waitbar(0,'最近邻域插值');
% for i = 1:Na
%     for j = 1:Nr
%         fai = fdoc + (i - Na/2) / Na * PRF;%方位向频率
%         d = sqrt(1 - (lambda*fai/2/v)^2);%D
%         r0 = (yc + (j - Nr/2)*Ts*C/2)*cos(theta);%最短斜距
%         ksrc = 2*v^2*fc^3/C/r0*d^3/fai^2;%Ksrc
%         
%         rcm = r0*(1/d - 1);%rcm值 对应距离
%         n_rcm = 2*rcm/C*Fs;%对应距离门
%         
%         delta_nrcm = n_rcm - floor(n_rcm);%小数部分
%         
%         if j + round(n_rcm) > Nr
%             signal_rcmc(i,j) = signal_RD(i,Nr/2);
%         else
%             if delta_nrcm >= 0.5
%                 signal_rcmc(i,j) = signal_RD(i,j+ceil(n_rcm));
%             else
%                 signal_rcmc(i,j) = signal_RD(i,j+floor(n_rcm));
%             end
%         end
%     end
%     waitbar(i/Na);
% end
% close(win);

win = waitbar(0,'sinc 8点插值');
corelen = 8;
signal_rcmc = zeros(Na,Nr);

for i = 1:Na
    for j = corelen:Nr
        fai = fdoc + (i - Na/2) / Na * PRF;%方位向频率
        d = sqrt(1 - (lambda*fai/2/v)^2);%D
        r0 = (yc + (j - Nr/2)*Ts*C/2)*cos(theta);%最短斜距
%         ksrc = 2*v^2*fc^3/C/r0*d^3/fai^2;%Ksrc
        
        rcm = r0*(1/d - 1);%rcm值 对应距离
        n_rcm = 2*rcm/C*Fs;%对应距离门
        
        delta_nrcm = n_rcm - floor(n_rcm);%小数部分
        
        for k = -corelen/2:corelen/2-1
            if n_rcm+k+j > Nr
                signal_rcmc(i,j) = signal_rcmc(i,j) + signal_RD(i,Nr) * sinc(k+n_rcm);
            else
                signal_rcmc(i,j) = signal_rcmc(i,j) + signal_RD(i,j+floor(n_rcm)+k) * sinc(k+delta_nrcm);
            end
        end
    end
       waitbar(i/Na);
end
close(win);

fu = linspace(fdoc - PRF/2,fdoc + PRF/2,Na);%方位向频率
% f = fftshift(fu);
f = fu;
d = sqrt(1 - (lambda*f/2/v).^2); 
ref_A = exp(1i*4*pi/lambda*(ones(Na,1)*range) .* (d'*ones(1,Nr))) ;%方位脉压滤波器
final_signal = iftx(signal_rcmc .* ref_A) .* (hamming(Na) * ones(1,Nr));

figure;
subplot(211);
imagesc(abs(echo));
xlabel('距离向');
ylabel('方位向');
title('回波信号');
subplot(212);
imagesc(abs(signal_comp));
xlabel('距离向');
ylabel('方位向');
title('距离脉压之后的信号');

figure;
subplot(211);
imagesc(abs(signal_RD));
xlabel('距离向');
ylabel('方位向');
title('未进行RCMC距离多普勒信号');
subplot(212);
imagesc(abs(signal_rcmc));
xlabel('距离向');
ylabel('方位向');
title('RCMC之后的距离多普勒信号');

figure;
imagesc(abs(final_signal));
xlabel('距离向');
ylabel('方位向');
title('最终的点目标');

figure;
mesh(abs(final_signal));
xlabel('距离向');
ylabel('方位向');
title('最终的点目标');

figure;
plot(20*log10(abs(final_signal) / abs(max(max(final_signal)))));
title('峰值旁瓣比');