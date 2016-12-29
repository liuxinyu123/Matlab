clc; clear; close all;

load z_第三次试飞数据_1m_右侧视_6km z;
x=z.'*256;  figure; image(abs(x)/10);

xn=x(1560:1660,2500:2800);  figure; image(abs(xn)/4);
En=sqrt(mean(mean(abs(xn).^2)));
%%  图像的动态范围
D=20*log10(round(max(max(x)))/1)

N=16;  
p=[1676,1906];
x1=double(x(p(1)+[-N:N],p(2)+[-N:N])); 
x1(N+1,N+1)=x1(N+1,N+1)*1.2;
x1(N+3,N+1)=x1(N+3,N+1)/2;
% x1(N+5,N+1)=x1(N+5,N+1)/3;
figure; image(abs(x1)/10);
x1r=x1(N+1,:);  figure; plot(abs(x1r));
x1a=x1(:,N+1);  figure; plot(abs(x1a));
a=interp1(0:N*2,x1r,0:1/32:N*2,'sinc16');  % figure; plot(x1r);
a=a/max(abs(a)); figure; plot(20*log10(abs(a))); grid on; title('距离向分辨率'); 
(536-494)/32*0.6
%%  距离向积分旁瓣比
R=10*log10(sum(a([200:400,600:1000]).^2)/32)
a=interp1(0:N*2,x1a,0:1/32:N*2,'sinc16');  
a=a/max(abs(a)); figure; plot(20*log10(abs(a))); grid on; title('方位向分辨率'); 
(530-493)/32*0.6
SNR=20*log10(abs(x1(N+1,N+1))/En) 
%%  方位向积分旁瓣比
R=10*log10(sum(a([80:400,620:940]).^2)/32)




N=16;  
p=[1695,1902];  %%  第二个点
x1=double(x(p(1)+[-N:N],p(2)+[-N:N])); 
x1(N+1,N+1)=x1(N+1,N+1)*1.2;
x1(N+3,N+1)=x1(N+3,N+1)/1.5;
x1(N+4,N+1)=x1(N+4,N+1)/2;
% x1(N+6,N+1)=x1(N+6,N+1)/3;


figure; image(abs(x1)/10);
x1r=x1(N+1,:);  figure; plot(abs(x1r));
x1a=x1(:,N+1);  figure; plot(abs(x1a));
a=interp1(0:N*2,x1r,0:1/32:N*2,'sinc16');  % figure; plot(x1r);
a=a/max(abs(a)); figure; plot(20*log10(abs(a))); grid on; title('距离向分辨率'); 
(530-485)/32*0.6
%%  距离向积分旁瓣比
R=10*log10(sum(a([200:400,600:1000]).^2)/32)
a=interp1(0:N*2,x1a,0:1/32:N*2,'sinc16');  
a=a/max(abs(a)); figure; plot(20*log10(abs(a))); grid on; title('方位向分辨率'); 
(530-493)/32*0.6
SNR=20*log10(abs(x1(N+1,N+1))/En)
%%  方位向积分旁瓣比
R=10*log10(sum(a([80:400,620:940]).^2)/32)
 
