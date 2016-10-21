%% Chirp Scaling算法
%徐一凡
clear all;clc;
%%距离向参数range:x domain
Tr=200;%时宽200m
Br=1;%带宽1
Kr=Br/Tr;%调频斜率
Fc=4;%载频4
Nfast=512;%为了快速运算
Xc=1200;X0=150;%定义距离向范围
x=Xc+linspace(-X0,X0,Nfast);%x域序列:Xc-X0~Xc+X0
dx=2*X0/Nfast;%定义步长
kx=linspace(-1/dx/2,1/dx/2,Nfast);%kx域序列

%%方位向参数cross-range:y domain
Ta=300;%时宽300m,合成孔径长度
Ba=1;%带宽1(1/m)
Ka=Fc/Xc;%调频斜率 Ka=Ba/Ta=Fc/Xc
Nslow=1024;%为了快速运算
Y0=200;
y=linspace(-Y0,Y0,Nslow);%y域序列:-Y0~Y0
dy=2*Y0/Nslow;
ky=linspace(-1/dy/2,1/dy/2,Nslow);%ky域序列

%%目标几何关系target geometry
%x坐标,y坐标,复后向散射系数 
Ptar=[Xc,0,1+0j              
          Xc+50,-50,1+0j
          Xc+50,50,1+0j
          Xc-50,-50,1+0j
          Xc-50,50,1+0j];
      
disp('Position of targets');disp(Ptar)
%%生成SAR正交解调后的回波数据
Srnm=zeros(Nfast,Nslow);
N=size(Ptar,1);%目标个数
h = waitbar(0,'SAR回波生成');
for i=1:1:N
    xn=Ptar(i,1);yn=Ptar(i,2);sigma=Ptar(i,3);%提取每个目标的信息
    X=x.'*ones(1,Nslow);%扩充为矩阵
    Y=ones(Nfast,1)*y;%扩充为矩阵
    DX=sqrt(xn^2+(Y-yn).^2);%中间变量
    phase=pi*Kr*(X-DX).^2-2*pi*Fc*DX;%回波相位
    Srnm=Srnm+sigma*exp(j*phase).*(abs(X-DX)<Tr/2).*(abs(Y-yn)<Ta/2);%回波累加
    waitbar(i/N)
end
close(h)
tic;
%%数据准备
phi0=-x'*sqrt(Fc^2-ky.^2);
phi1=-Fc*x'*(1./sqrt(Fc^2-ky.^2));
phi2=1/2*x'*(ky.^2./(Fc^2-ky.^2).^1.5);
Cs=ones(Nfast,1)*(Fc./sqrt(Fc^2-ky.^2)-1);
Ks=1./(1/Kr-2*phi2);

%%CSA:7步  开始
s_xky=fftshift(fft(fftshift(Srnm).')).';%方位向FFT(步骤1)
scs_xky=s_xky.*exp(j*pi*Cs.*Ks.*(x'*ones(1,Nslow)-Xc*(1+Cs)).^2);%Chirp Scaling（步骤2）
s1=ifty(scs_xky);%为显示存储数据
scs_kxky=fftshift(fft(fftshift(scs_xky)));%距离向FFT（步骤3）
srmc_kxky=scs_kxky.*exp(j*pi*(kx.^2'*ones(1,Nslow))./(1+Cs)./Ks...
                     +j*2*pi*Xc*Cs.*(kx'*ones(1,Nslow)));%距离迁移校正&距离向匹配滤波（步骤4）
srmc_xky=fftshift(ifft(fftshift(srmc_kxky)));%距离向IFFT（步骤5）
f_xky=srmc_xky.*exp(-j*pi*Ks.*Cs.*(1+Cs).*((x-Xc).^2'*ones(1,Nslow))...
           -j*2*pi*phi0);%消除误差函数，方位向匹配滤波（步骤6）
f_xy=fftshift(ifft(fftshift(f_xky).')).';%方位向IFFT（步骤7）
%%CSA:7步  结束
toc;

%%为了显示CSA算法的一致RCMC，将s1进行距离向压缩显示
p0_x=exp(j*pi*Kr*(x-Xc).^2).*(abs(x-Xc)<Tr/2);%距离向LFM信号
p0_kx=fftshift(fft(fftshift(p0_x)));
p0_y=exp(-j*pi*Ka*y.^2).*(abs(y)<Ta/2);%方位向LFM信号
p0_ky=fftshift(fft(fftshift(p0_y)));
 
s_kxy=fftshift(fft(fftshift(s1)));%距离向FFT
sxc_kxy=s_kxy.*(conj(p0_kx).'*ones(1,Nslow));
sxc_kxky=fftshift(fft(fftshift(sxc_kxy).')).';%距离压缩后的2D频域信号
sxc_xy=fftshift(ifft(fftshift(sxc_kxy)));%距离压缩后的信号
sxc_xky=fftshift(fft(fftshift(sxc_xy).')).';%距离压缩后，距离-多普勒域
%%结果显示
figure(1)
colormap(gray);
imagesc(255-abs(Srnm));
xlabel('方位向'),ylabel('距离向'),
title('仿真出来的信号');

figure(2)
colormap(gray);
imagesc(255-abs(sxc_xy));
xlabel('方位向'),ylabel('距离向'),
title('Chirp Scaling后、经过距离向压缩，距离徙动一致');

figure(3)
colormap(gray);
imagesc(255-abs(srmc_xky)); 
xlabel('方位向'),ylabel('距离向'),
title('消除距离徙动后的信号');

figure(4)
colormap(gray);
imagesc(255-abs(f_xky)); 
xlabel('方位向'),ylabel('距离向'),
title('相位校正后的信号');

figure(5)
colormap(gray);
imagesc(255-abs(f_xy)); 
xlabel('方位向'),ylabel('距离向'),
title('生成的点目标');
