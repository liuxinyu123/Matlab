%%CS 算法 by W&Z---------------------------------已看
%%================================================================
clear;clc;close all;
%%================================================================
%%Parameter--constant
C=3e8;                                   %propagation speed
%%Parameter--radar characteristics
Fc=10e9;                                  %carrier frequency 10GHz 发射机频率
lambda=C/Fc;                             %wavelength                
%%Parameter--target area
Xmin=-50;                                  %target area in azimuth is within[Xmin,Xmax]方位向
Xmax=50;                                  
Yc=10000;                                %center of imaged area
Y0=500;                                  %target area in range is within[Yc-Y0,Yc+Y0]距离向
                                         %imaged width 2*Y0
%%Parameter--orbital information
V=100;                                   %SAR velosity 100 m/s 正侧视
H=5000;                                  %height 5000 m
R0=sqrt(Yc^2+H^2);                       %目标离飞机的最近距离        
%%Parameter--antenna
D=4;                                     %antenna length in azimuth direction方位向天线长度
Lsar=lambda*R0/D;                        %合成孔径长度    
Tsar=Lsar/V;                             %合成孔径时间         

%% 方位向
Ka=-2*V^2/lambda/R0;                      %方位向多普勒斜率
Ba=abs(Ka*Tsar);                          %方位向多普勒调频范围即多普勒带宽
PRF=Ba;                                   %方位向发射序列频率就等于方位向多普勒频带宽度？？？？？从何而来
PRT=1/PRF;                                %计算方位向发射序列间隔
ds=PRT;                                   %脉冲重复时间
Nslow=ceil((Xmax-Xmin+Lsar)/V/ds);        %慢时间采样点数，但是考虑到边缘，须多加一个合成孔径的采样时间（方位向变换时间/脉冲重复时间）
Nslow=2^nextpow2(Nslow);                  %for fft
sn=linspace((Xmin-Lsar/2)/V,(Xmax+Lsar/2)/V,Nslow);     %产生方位向信号离散序列，以时间为量纲
PRT=(Xmax-Xmin+Lsar)/V/Nslow;            %方位向单位间隔
PRF=1/PRT;                               %采样频率
ds=PRT;                                  %采样间隔（时间域） 

%% 距离向
Tr=5e-6;                                 %pulse duration 5us
Br=30e6;                                 %chirp frequency modulation bandwidth 30MHz
Kr=Br/Tr;                                %调频斜率，正扫频
Fsr=3*Br;                                %采样频率
dt=1/Fsr;                                %采样间隔
Rmin=sqrt((Yc-Y0)^2+H^2);                %最小斜距
Rmax=sqrt((Yc+Y0)^2+H^2+(Lsar/2)^2);     %最大斜距                                  
Nfast=ceil(2*(Rmax-Rmin)/C/dt+Tr/dt);    %快时间采样点数
Nfast=2^nextpow2(Nfast);                 %for fft
tm=linspace(2*Rmin/C,2*Rmax/C+Tr,Nfast); %产生距离向信号离散序列，以时间为量纲                     
dt=(2*Rmax/C+Tr-2*Rmin/C)/Nfast;         %距离向采样间隔
Fsr=1/dt;                                %距离向采样频率 
%%
%%Parameter--resolution                      
DY=C/2/Br;                                %距离向分辨率  
DX=D/2;                                   %方位向分辨率
%%Parameter--point targets
Ntarget=3;                                %三个点目标
%format [x, y, reflectivity]
Ptarget=[Xmin+10*DX,Yc+10*DY,1               
              Xmin+30*DX,Yc,1
              Xmin+20*DX,Yc+20*DY,1];     %目标和参考点的位置矩阵及其散射系数
%%================================================================
%%Generate the raw signal data
K=Ntarget;                                 %number of targets
N=Nslow;                                   %number of vector in slow-time domain  %方位向 
M=Nfast;                                   %number of vector in fast-time domain  %距离向
T=Ptarget;                                 %position of targets
Srnm=zeros(N,M);
for k=1:1:K
    sigma=T(k,3);                           %目标散射系数
    Dslow=sn*V-T(k,1);                      %各个方位向点到目标的方位向距离
    R=sqrt(Dslow.^2+T(k,2)^2+H^2);          %目标离各个方位向上的点的斜距   
    tau=2*R/C;                              %回波延时
    Dfast=ones(N,1)*tm-tau'*ones(1,M);      %线性调频分量               
    phase=pi*Kr*Dfast.^2-(4*pi/lambda)*(R'*ones(1,M));%回波相位=它等于多普勒分量加上线性调频分量，并且要加窗处理？？？？？？？？？？？？相位符号？
    Srnm=Srnm+sigma*exp(1i*phase).*(0<Dfast&Dfast<Tr).*((abs(Dslow)<Lsar/2)'*ones(1,M));
end
%%================================================================


%%
%% CS参数设置
f=linspace(-PRF/2,+PRF/2,Nslow);             %方位向
fr=linspace(-Fsr/2,+Fsr/2,Nfast);            %距离向
r=tm/2*C;                                    %距离
r_ref=sqrt(Yc^2+H^2);
CS_f=1./sqrt(1-(lambda*f/2/V).^2)-1;    %%%---------------（7.17式）1/D-1
% Ks=Kr./( 1+Kr*r_ref*(2*lambda/C^2)*((lambda*f/2/V).^2)./(1-(lambda*f/2/V).^2).^1.5 ); %%--------（7.18式）
Ks=Kr./( 1-Kr*r_ref*(2*lambda/C^2)*((lambda*f/2/V).^2)./(1-(lambda*f/2/V).^2).^1.5 ); %%--------（7.18式）
R_ref=r_ref*(1+CS_f);                   %%------------ Rref/D
%% CS成像处理，chirp scaling变换
Srnm_xfft=fftshift(fft(fftshift(Srnm)));
% CS_phase=exp( -j*pi*(Ks.*CS_f)'*ones(1,M).*( ones(N,1)*tm-(2*R_ref/C)'*ones(1,M) ).^2 );%%----------（7.30式）
CS_phase=exp( 1i*pi*(Ks.*CS_f)'*ones(1,M).*( ones(N,1)*tm-(2*R_ref/C)'*ones(1,M) ).^2 );%%----------（7.30式）
Srnm_cs=Srnm_xfft.*CS_phase;                   %CS相位因子相乘

figure;
mesh(abs(Srnm_cs));
%% 一致RCMC和距离压缩

Srnm_yfft=fftshift(fft(fftshift(Srnm_cs.'))).';
% rmc_phase=exp( -j*pi*(ones(N,1)*fr.^2)./((Ks.*(1+CS_f))'*ones(1,M)) ).*exp( j*4*pi/C*(ones(N,1)*fr)*r_ref.*(CS_f'*ones(1,M)) );%%------（7.32式）
rmc_phase=exp(1i*pi*(ones(N,1)*fr.^2)./((Ks.*(1+CS_f))'*ones(1,M)) ).*exp( 1i*4*pi/C*(ones(N,1)*fr)*r_ref.*(CS_f'*ones(1,M)) );%%------（7.32式）
Srnm_rmc=Srnm_yfft.*rmc_phase;                 %距离徙动校正与距离压缩

figure;
mesh(abs(Srnm_rmc));

Srnm_yifft=fftshift(ifft(fftshift(Srnm_rmc.'))).';
r_sub=ones(N,1)*(r-r_ref).^2;
phase_cor=(4*pi/C^2*Ks.*(1+CS_f).*CS_f)'*ones(1,M).*r_sub;
% Srnm_xphase=exp( j*4*pi/lambda*ones(N,1)*r.*( sqrt(1-((lambda*f/2/V).^2)'*ones(1,M)) )-j*phase_cor-j*2*pi*C*ones(N,1)*tm/lambda   ); %%-----（7.34式）        
Srnm_xphase=exp(1i*4*pi/lambda*ones(N,1)*r.*( sqrt(1-((lambda*f/2/V).^2)'*ones(1,M)) )-1i*phase_cor-1i*2*pi*C*ones(N,1)*tm/lambda   ); %%-----（7.34式） 
% Srnm_xphase=exp(1i*4*pi/lambda*ones(N,1)*r.*( sqrt(1-((lambda*f/2/V).^2)'*ones(1,M)) )+1i*phase_cor);
Srnm_cor=Srnm_yifft.*Srnm_xphase;             %方位压缩与相位补偿

figure;
mesh(abs(Srnm_cor));
f_xy=fftshift(ifft(fftshift(Srnm_cor))); 
%% =======================================================================
%%plot image
figure;
Ga=abs(f_xy);
%figure,mesh(Ga((200:300),(750:860)));axis tight;
mesh(Ga);
figure
imagesc(Ga)






