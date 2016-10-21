nmb = 4;                               %目标个数4个   
rrec =200;                               %测量目标最远距离
b = 100e6;                               %调频信号带宽
smb_range = [10,30,40,100];                %三点目标的距离     最小分辨距离为s=c/2b=1.5m
smb_rcs = [1 1 1 2];                      %三点目标的横截面积
taup = 5e-6;                         %信号持续脉宽
f0 = 5e9;                              % 载频频率
c = 3e8;                                % 信号传播的速度，即光速
fs = 2*b;                                % 采样的频率
sampling_interval = 1/fs;                
n = fix(taup/sampling_interval);          %总共点数（取整）
nfft =n;                                   % 采样点数
freqlimit = 0.5*fs;                       
freq = linspace(-freqlimit,freqlimit,n);  % 频率采样间隔 = fs/n = 1/taup;
t = linspace(-taup/2,taup/2,n);            %相邻点时间间隔
x(:,1:n) = 0.;                             % x为矩阵
y(1:n) = 0.;     
replica(1:n) = 0.;
replica = exp(1i * pi * (b/taup) .* t.^2);  %基带线性调频信号
for j = 1:1:nmb                            %矩阵方法将接收信号叠加
    range = smb_range(j) ;
    x(j,:) = smb_rcs(j) .*exp(-1i*2*pi*f0*2*range/c).* exp(1i * pi * (b/taup) .* (t + 2*range/c).^2) ; %接收信号
    y = x(j,:)  + y;                       %信号叠加
end
rfft = fft(replica,nfft);
yfft = fft(y,nfft);
out= abs(ifft((rfft .* conj(yfft)))) ./ (nfft);
s = taup * c /2; Npoints = ceil(rrec * nfft /s);
dist =linspace(0, rrec, Npoints);          
%图片显示：
figure 
subplot(311)
plot(t,real(replica));axis tight;
xlabel('Range in meters');ylabel('Amplitude in dB');
title('线性调频信号');
subplot(312)
plot(t,real(y));axis tight;
xlabel('Range in meters');ylabel('Amplitude in dB');
title('压缩前雷达回波');
subplot(313);
plot(dist, out(1:Npoints));
xlabel ('Target relative position in meters');
ylabel ('压缩后雷达回波');
grid on; 