clear;
clc;
close all;
C=3e8;                                  %光速  
fc=5.3e9;                              %载频
lambda=C/fc;                            %波长
Da=3.75;                                   %天线孔径
theta_bw = lambda/Da;                   %方位波束宽度
Rs=20e3;                                %中心斜距(单位：m)  
squint_angle=0/180*pi;                  %天线波束斜视角，正侧视情况下为0度 
beta=60/180*pi;                         %下视角60度
H = Rs*cos(beta);                       %载机高度
Yb= Rs*sin(beta);
V=150;                                  %(机载)雷达有效速度(单位：m/s)
Ta=Rs*theta_bw/V;                       %合成孔径时间
Tp=10e-6;                               %发射脉冲时宽（单位：s） 
B=40e6;                                 %带宽
Fs=60e6;                               %距离采样率（单位：Hz)
gama=B/Tp;                             %距离向调频率（单位：Hz/s)
f_doppler = 2*V/lambda;                    %多普勒带宽
PRF=100;                               %方位采样率（单位：Hz)
nrn=floor(1.2*Tp*Fs/2)*2;          %距离采样点数
nan=floor(1.2*Ta*PRF/2)*2;                 %距离线数（也就是方位采样点数）
nrn=1024;
nan=200;
x=zeros(nrn,nan);                       %定义雷达回波数据矩阵
tnrn=[-nrn/2:nrn/2-1]'/Fs+2*Rs/C;     %距离向时间数据，加上2Rs/c是为了定位景中心...
                                      ...目标,由c*t/2=Rs得来  ？？？？        加上2Rs/c  目的
%tnrn=[-nrn/2:nrn/2-1]'/Fs; 
Ryta=tnrn*C/2;
tnan=[-nan/2:nan/2-1]'/PRF;              %方位向时间数据
PointN=1;
sigman=1;
x0=0;y0=Yb;z0=0;
theta_Geo=theta_bw;
P=zeros(nrn,nan);
for n=1:nan
    xT=V*tnan(n);
    yT=0;
    zT=H;
    xn=x0-xT;                            %成像点位置与雷达位置距离的x方向投影
    yn=y0-yT;                            %成像点位置与雷达位置距离的y方向投影
    zn=z0-zT;
    R=sqrt(xn^2+yn^2+zn^2);            %成像点位置与雷达位置距离的大小
    Ka=2*V^2/lambda./Ryta;
        if abs(atan(xn/yn))<theta_Geo/2        %判断点目标是否在雷达的波束角范围内
        winr=((tnrn-2*R/C)<=Tp/2)&((tnrn-2*R/C)>=-Tp/2);    %回波是否在发射信号包络中       ？？？？？？？    1/4、3/4分割行不行？？
        echo=winr.*exp(j*pi*gama*(tnrn-2*R/C).^2)*exp(-j*4*pi*R/lambda);
        x(:,n)=x(:,n)+echo;
        P(:,n)=winr.*(-pi*Ka*tnan(n)^2+pi*gama*(tnrn-2*R/C).^2);
        %P(:,n)=winr.*(-pi*Ka*tnan(n)^2+pi*gama*(tnrn-2*R/C).^2);
        end

        
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%原始信号%%%%%%%%%%%%%%%%%%%%%%%%%
%P=(abs(x)>0).*P;
figure,subplot(221),imagesc(abs(x));
title('(a)回波信号幅度');
xlabel('方位向');                 %这一部分是画原始回波数据的时域幅度图
ylabel('距离向');
%subplot(222),contour(P,60);
title('(b)回波信号相位');
xlabel('方位向');                 %这一部分是画原始回波数据的时域幅度图
ylabel('距离向');
save x;
%wr=(kaiser(nrn,2.5))';                    %对匹配滤波器加Kaiser窗，可根据需要加...
                                       ...Taylor/Chebyshev/Hanning/Hamming/Kaiser窗等
tp_nrn=floor(Tp*Fs/2)*2;
hrt=exp(-j*pi*gama*([-tp_nrn/2:tp_nrn/2-1]'/Fs).^2);
serf=zeros(nrn,1);
serf(nrn/2-tp_nrn/2:nrn/2+tp_nrn/2-1,1)=hrt;
Hrf=fft(serf);                            %把匹配滤波器变换到频域
for n=1:nan
x(:,n)=fftshift(ifft(fft(x(:,n)).*Hrf));           %这一部分实现的是距离向的匹配滤波，即距离压缩
end 
subplot(223),imagesc(abs(x));
title('距离压缩');
xlabel('方位向(采样点)');                 %这一部分是画距离压缩后的幅度图
ylabel('距离向(采样点)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%方位压缩
for m=1:nrn
   % tempFlag = 0;
   % //if tempFlag==1
    %Ka=2*V^2/lambda/Ryta(m);        %多普勒调制频率    ???????????   为什么不是2式,为什么不是负的
    %else
    Ka=2*V^2/lambda/Rs;        %多普勒调制频率 2式
    %end
    hat=(exp(-j*pi*Ka*tnan.^2))';     %???  
    Haf=fft(hat);
    x(m,:)=fftshift(ifft(fft(x(m,:)).*Haf));
end
subplot(224),imagesc(abs(x));
title('方位压缩');
xlabel('方位向(采样点)');                 %这一部分是画方位压缩后的幅度图
ylabel('距离向(采样点)');