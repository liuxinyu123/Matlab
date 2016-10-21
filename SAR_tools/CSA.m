
clear all;   close all;
t1=cputime;                                                  
Fs=32.317e+006;                            
Fr=7.2135e+11   ;                     
start=6.5956e-003;                     
Tr=4.175e-05 ;                               
R0=9.88646462e+05   ;                        
f0=5.3e+09  ;                          
lamda=0.05667;                    
Fa=1.25698e+03 ;                    
Vr=7062;                             
Kr=0.72135e+012;
Ka=1733;                        
Fc=-6900;                             
c=299790000;                    
load CDdata1.mat
data=double(data);                        
[length_a,length_r]=size(data);                
T_start=6.5956e-003;                   %数据窗开始时间
tau=T_start:double(1/Fs):T_start+double(length_r/Fs)-double(1/Fs);  
R_ref=(T_start+length_r/Fs)/2*c;            
f_a=(-Fa/2+Fc):(Fa/length_a):(Fa/2+Fc-Fa/length_a);      
f_r=0:Fs/length_r:Fs-Fs/length_r;            
D = (1 - (f_a*lamda/2/Vr).^2).^0.5;         %距离徙动因子
alpha = 1./D - 1;                    
R = R_ref./D;                           %距离多普勒域中更精确的双曲线距离等式
Z=(R0*c*f_a.^2)./(2*Vr^2*f0^3.*D.^3);
Km=Kr./(1-Kr.*Z);                             %校正后距离脉冲调频率
for i=1:length_r
    data(:,i)=fft(data(:,i));
end
tau1=ones(length_a,1)*tau;
tau2=2.*(R(:)*ones(1,length_r))./c;
D_tau=tau1-tau2;
H1=Km.*alpha;
H1=H1(:)*ones(1,length_r);
Ssc=exp(-j*pi*H1.*D_tau.^2);               % 线性变标方程
data=data.*Ssc;                           %校正补余RCM
for i=1:length_a
    data(i,:)=fftshift(fft(fftshift(data(i,:))));
end
H_r_1=1./(Km.*(1+alpha));
H_r_1=H_r_1(:)*ones(1,length_r);
H_r_2=ones(length_a,1)*f_r.^2;
H_r=exp(-j*pi*H_r_1.*H_r_2);
H_RCMC=exp(j*4*pi*R_ref.*(ones(length_a,1)*f_r).*(alpha(:)*ones(1,length_r))/c); 
data=data.*H_r.*H_RCMC;                       %距离压缩和一致RCMC
for i=1:length_a
    data(i,:)=fftshift(ifft(fftshift(data(i,:))));
end
r_0=tau/2*c;
H3_1=Km.*alpha.*(1+alpha)./(c^2);
H3=H3_1(:)*ones(1,length_r);
phi=4*pi.*H3.*(ones(length_a,1)*(r_0-R_ref).^2)/(c^2);
data=data.*exp(j*phi);                                % 完成相位校正
H_a= exp(j*4*pi/lamda*(ones(length_a,1)*r_0).*((D(:)-1)*ones(1,length_r)));
data = data .*H_a;                 % 方位向匹配滤波
for i=1:length_r
    data(:,i)=ifft(data(:,i));
end
data = fftshift(data,2);           
figure;
colormap(gray(256))
imagesc(log(10*abs(data)));xlabel('Range');ylabel('Azimuth');
Compute_time=cputime-t1               %输出运行时间

