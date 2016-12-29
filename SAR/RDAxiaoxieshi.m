%2016/11/8
%Liu Yakun

clc,clear,close all;
C = 3e8;%光速
lambda = 0.03;            %发射信号波长
fc = C / lambda;          %发射信号中心频率
v = 200;%雷达平台速度
% h = 5000;
D = 5;%天线长度
theta = 15 / 180 * pi;%斜视角
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
 
%%%%%%%%%%%%%%%%%%%Range Compression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_range = [-Nr/2:Nr/2-1]/Nr*Fs;
H_range = exp(1i*pi*f_range.^2/kr) .* kaiser(Nr,2.5)';
signal_comp = ifty(fty(echo) .* (ones(Na,1) * H_range));
