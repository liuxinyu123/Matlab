%2016/12/1
%Liu Yakun

clc,clear,close all;
%%%%%%%%%%%%%%%%%%%%%%参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = 3e8;                  %光速
Theta = 45 / 180 * pi;    %斜视角
Rc = 41.7e3;              %场景中心斜距
R0 = Rc * cos(Theta);     %场景中心最近距离
lambda = 0.03;            %发射信号波长
fc = C / lambda;          %发射信号中心频率
Tr = 2e-6;                %脉冲 持续时间
Br = 60e6;                %带宽
Kr = Br / Tr;             %距离调频率
V = 200;                  %平台飞行速度
D = 5;                    %方位向天线长度
Beta = lambda / D;        %波束宽度
Lsar = lambda * Rc / D;   %正侧视合成孔径长度
Lsar_squint = R0 *(tan(Theta + Beta/2) - tan(Theta - Beta/2));
Tsar = Lsar_squint / V;
Dx = D / 2;
Dy = C / 2 / Br;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%场景参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scene_center = [Rc*cos(Theta) Rc*sin(Theta)];
range_width = 500;
azimuth_width = 500;
Ymin = scene_center(1) - range_width/2;
Ymax = scene_center(1) + range_width/2;
Xmin = scene_center(2) - azimuth_width/2;
Xmax = scene_center(2) + azimuth_width/2;


targets = [ scene_center(1)       scene_center(2)        1
%             scene_center(1)+200   scene_center(2)        1
%             scene_center(1)-200   scene_center(2)        1
%             scene_center(1)       scene_center(2)+200  1
%             scene_center(1)       scene_center(2)-200  1
%             scene_center(1)+200   scene_center(2)+200  0.7
%             scene_center(1)-200   scene_center(2)-200  1
%             scene_center(1)+200   scene_center(2)-200  0.9
%             scene_center(1)-200   scene_center(2)+200  0.4
            ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%显示场景中的点坐标%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(targets(:,1),targets(:,2),'b*');
xlabel('Range axis [m]');
ylabel('Azimuth axis [m]');
title('Targets in scene');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%慢时间%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ta_start = (Xmin - Ymax * tan(Theta + Beta/2)) / V;
ta_end = (Xmax - Ymin * tan(Theta - Beta/2)) / V;

Ka = -2*V^2*cos(Theta)^2/lambda/Rc;
Ba = abs(Ka*Tsar);
alpha_PRF = 1.3;

PRF =  2 * round(alpha_PRF * Ba);
Na = (ta_end - ta_start) * PRF; 
Na = 2^nextpow2(Na);

fdoc = 2 * V * sin(Theta) / lambda;
% Mamb = floor(fdoc / PRF);
% fdoc = fdoc - Mamb*PRF;
Tslow = [-Na/2:Na/2 - 1] / PRF;
distance_azimuth = Tslow * V;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%快时间%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rnear = Ymin / cos(Theta - Beta/2);
Rfar = Ymax / cos(Theta + Beta/2);

PRFmin = Ba + 2*V*Kr*Tr*sin(Theta)/C;
RFmax = 1/(2*Tr+2*(Rfar-Rnear)/C);
Rmid = (Rnear + Rfar) / 2;
alpha_Fr = 1.2;
Fr = round(alpha_Fr * Br);
Nr = (2*(Rfar - Rnear)/C + Tr) * Fr;
Nr = 2^nextpow2(Nr);

Tfast = [-Nr/2:Nr/2 - 1]/Fr + 2 * Rmid / C;
distance_range = Tfast * C / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%回波产生%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo = zeros(Na,Nr);
nTar = size(targets,1);
disp('*****************************************');
disp('参数打印:');
disp(['慢时间采样点数为:',poly2str(Na),'点']);
disp(['脉冲重复发射频率为:',poly2str(PRF),'Hz']);
disp(['快时间采样点数为:',poly2str(Nr),'点']);
disp(['快时间采样频率为:',poly2str(Fr),'Hz']);
disp(['方位向合成孔径长度为:',poly2str(Lsar),'m']);
disp(['方位向分辨率为:',poly2str(Dx),'m']);
disp(['距离向分辨率为:',poly2str(Dy),'m']);
disp('*****************************************');
disp('产生回波信号.....');
for i = 1:nTar
    range = sqrt(targets(i,1)^2 + (targets(i,2) - V*Tslow).^2);
    tau = 2 * range / C;
    D = ones(Na,1) * Tfast - tau.' * ones(1,Nr);
    phase = pi*Kr*D.^2 - 4*pi/lambda*(range.' * ones(1,Nr));
    sigma = targets(i,3);
    tci = (targets(i,2) - scene_center(2) + tan(Theta) * (targets(i,1) - scene_center(1))) / V;
    tsi = targets(i,1) * (tan(Theta + Beta/2) - tan(Theta)) / V;
    tei = targets(i,1) * (tan(Theta) - tan(Theta - Beta/2)) / V;
%     time_wave_center = (targets(i,2) - scene_center(2))/V;
    echo = echo + sigma * exp(1i * phase) .* (abs(D) < Tr/2) .* (((Tslow > (tci - tsi) & Tslow < (tci + tei)))' * ones(1,Nr));
end

figure;
imagesc(abs(echo));
xlabel('Range axis');
ylabel('Azimuth axis');
title('Raw signal');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%Range Chirp signal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% chirp_range = zeros(1,Nr);
% t_range = -Tr/2:1/Fr:Tr/2;
% omega_range = -Fr/2:1/Tr:Fr/2;
% chirp_range_temp = exp(1i * pi * Kr * t_range.^2);
% length_range = length(t_range);
% chirp_range(ceil((Nr - length_range) / 2):ceil((Nr - length_range) / 2) + length_range - 1) = chirp_range_temp;
% 
% figure;
% subplot(211);
% plot(t_range,real(chirp_range_temp),'r');
% hold on;
% plot(t_range,imag(chirp_range_temp),'g');
% xlabel('Time axis [sec]');
% ylabel('Amplitude');
% xlim([min(t_range) max(t_range)]);
% title('Range Chirp');
% 
% subplot(212);
% plot(omega_range,abs(fftshift(fft(chirp_range_temp))));
% xlabel('Frequence [Hz]');
% ylabel('Absolute');
% xlim([min(omega_range) max(omega_range)]);
% title('Spectrum of range chirp');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Azimuth Chirp signal%%%%%%%%%%%%%%%%%%%%%%
% chirp_azimuth = zeros(1,Na);
% t_azimuth = -Tsar/2:1/PRF:Tsar/2;
% omega_azimuth = (fdoc - PRF/2):1/Tsar:(fdoc + PRF/2);
% chirp_azimuth_temp = exp(1i * pi * Ka * t_azimuth.^2);
% length_azimuth = length(t_azimuth);
% chirp_azimuth(ceil((Nr - length_azimuth) / 2):ceil((Nr - length_azimuth) / 2) + length_azimuth - 1) = chirp_azimuth_temp;
% 
% figure;
% subplot(211);
% plot(t_azimuth,real(chirp_azimuth_temp),'r');
% hold on;
% plot(t_azimuth,imag(chirp_azimuth_temp),'g');
% xlabel('Time axis [sec]');
% ylabel('Amplitude');
% xlim([min(t_azimuth) max(t_azimuth)]);
% title('Azimuth Chirp');
% 
% subplot(212);
% plot(omega_azimuth,abs(fftshift(fft(chirp_azimuth_temp))));
% xlabel('Frequence [Hz]');
% ylabel('Absolute');
% xlim([min(omega_azimuth) max(omega_azimuth)]);
% title('Spectrum of range chirp');
%%%%%%%%%%%%%%%%%%%%%%距离压缩%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t = Tfast - 2 * Rc / C;
% chirp_range = exp(1i * pi * Kr * t.^2) .* (abs(t) < Tr / 2);
% signal_Ra = fty(echo) .* (ones(Na,1) * conj(fty(chirp_range)));
disp('距离压缩开始.......');
f_range = [-Nr/2:Nr/2-1]/Nr*Fr;
f_azimuth = fdoc + [-Na/2:Na/2-1]/Na*PRF;
d = sqrt(1-(lambda*f_azimuth/2/V).^2);
Mamb = round(fdoc / PRF);
delta_fdoc = fdoc - Mamb * PRF;
% R0 = Rmid * cos(Theta);
% H = exp(-1i*4*pi/C*(ones(Na,1) * f_range) .* ((R0./d - R0)' * ones(1,Nr)));
% echo = echo .* H;

H_range = exp(1i * pi * (f_range.^2) / Kr) .* kaiser(Nr,2.5)';
alpha = [-Fr/2:Fr/2-1]/Fr;
signal_Ra = fty(echo) .* (ones(Na,1) * H_range) .* (exp(-1i*2*pi*delta_fdoc.*Tslow).' * ones(1,Nr));
delta_f = linspace(-20/2,20/2,Nr);
% for i = 1:Nr
%     signal_Ra(:,i) = signal_Ra(:,i) .* exp(-1i*2*pi*delta_f(i)*Tslow');
% end
signal_comp = ifty(signal_Ra);
disp('距离压缩结束');
% signal_rA = ftx(signal_comp);
% 
% figure;
% plot(abs(signal_comp1(512,:)),'r');
% hold on;
% plot(abs(signal_comp2(512,:)),'g');
%%%%%%%%%%%%%%%%%%%%%%%%%%显示距离压缩后的图像%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
colormap('gray');
imagesc(255 - abs(signal_comp));
xlabel('Range axis');
ylabel('Azimuth axis');
title('Range compression');

% signal_rA = ftx(signal_comp);
% delta_range = C/2/Fr; 
% R0 = Rmid * cos(Theta);
% n_rcm=(1./sqrt(1-(f_azimuth*lambda/2/V).^2)'-1)*R0/delta_range;
% 
% alpha_Nr=([-Nr/2:Nr/2-1]/Nr);
% for i = 1:Na
%     signal_rA(i,:) = fftshift(ifft(fftshift(fft(fftshift(signal_rA(i,:))) .* exp(1i*2*pi*n_rcm(i)*alpha_Nr))));
% end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SRC in RFAF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
% signal_RA = fty(signal_rA);
disp('SRC开始........');
signal_RA = ftx(signal_Ra);
r = R0;
h_src = exp(-1i*pi*r*C/2/V^2/fc^3*(f_azimuth'.^2 ./ d'.^3  * ones(1,Nr)) .* (ones(Na,1) * f_range.^2));
signal_RA_src = signal_RA .* h_src;
% R0 = Rmid * cos(Theta);
% H = exp(-1i*4*pi/C*(ones(Na,1) * f_range) .* ((R0./d - R0)' * ones(1,Nr)));
% signal_RA_src = signal_RA_src .* H;
signal_rA_src = ifty(signal_RA_src);
% signal_rA_src_shift = zeros(Na,Nr);
% signal_rA_src_shift = signal_rA_src()
% signal_Ra_src = iftx(signal_RA_src);
signal_ra_src = iftx(signal_rA_src);
disp('SRC结束');
% Mamb = round(fdoc / PRF);
% delta_fdoc = fdoc - Mamb * PRF;
% signal_rA_src = ftx(signal_ra_src .* (exp(-1i*2*pi*delta_fdoc*Tslow).' * ones(1,Nr)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RCMC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% range = distance_range * cos(Theta);
% delta_range = C/2/Fr;
% rcm = (ones(Na,1) * range) .* ((1./d' - 1) * ones(1,Nr));
% h_rcmc = exp(1i*4*pi/C*(ones(Na,1) *f_range) .* rcm);
% signal_RA_src_rcm = signal_RA_src .* h_rcmc;
% signal_rA_src_rcm = ifty(signal_RA_src_rcm);
% FF = ones(Na,1) * f_range;
% FU = f_azimuth' * ones(1,Nr);
% Refcorr = exp(1i * pi / fc^2 / Ka * (FU.*FF).^2 + 1i * pi * fdoc^2 / fc / Ka * FF - 1i * pi / fc / Ka * FU.^2 .* FF); 
% signal_RA_src_rcmc = signal_RA_src .* Refcorr;
% signal_rA_src_rcmc = ifty(signal_RA_src_rcmc);
delta_range = C/2/Fr; 
R0 = Rmid * cos(Theta);
n_rcm=(1./sqrt(1-(f_azimuth*lambda/2/V).^2)'-1)*R0/delta_range;
base_n_rcm = (1./sqrt(1-(f_azimuth(Na/2)*lambda/2/V).^2)-1)*R0/delta_range;
n_rcm = n_rcm - base_n_rcm;
signal_rA_src_rcmc = zeros(Na,Nr);
alpha_Nr=([-Nr/2:Nr/2-1]/Nr);
length = 8;
disp('距离徙动校正开始.......');
win = waitbar(0,'Sinc 8点插值');
for i = 1:Na
    for j = length:Nr
        rcm = n_rcm(i);
        if rcm > 0
            delta_n_rcm = rcm - floor(rcm);
        else            
            delta_n_rcm = rcm - ceil(rcm);
        end
        
        for k = -length/2:length/2 - 1
            if (k+j+ceil(rcm) > Nr) || (k+j+floor(rcm) <= 0)
                signal_rA_src_rcmc(i,j) = signal_rA_src_rcmc(i,j) + signal_rA_src(i,j) * sinc(k+rcm);
            else
                signal_rA_src_rcmc(i,j) = signal_rA_src_rcmc(i,j) + signal_rA_src(i,j+floor(rcm)+k) * sinc(delta_n_rcm+k);        
            end
        end
    end
    waitbar(i/Na);
end
close(win);
disp('距离徙动校正结束');
% for i = 1:Na
%     signal_rA_src(i,:) = fftshift(ifft(fftshift(fft(fftshift(signal_rA_src(i,:))) .* exp(1i*2*pi*n_rcm(i)*alpha_Nr))));
% end 

% figure;
% plot(abs((signal_rA_src_norcmc(512,:))),'r');
% hold on;
% plot(abs(signal_rA_src(512,:)),'g');

%%%%%%%%%%%%%%%%%%%%%Doppler Phase Compensation%%%%%%%%%%%%%%%%%%%%%%%%%
% x0 = linspace(-azimuth_width/2,azimuth_width/2,Nr) * tan(Theta);
% H_phase = exp(1i*2*pi/V*(f_azimuth .* x0));
% signal_rA_src = signal_rA_src .* (H_phase.' * ones(1,Nr));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Azimuth Compression%%%%%%%%%%%%%%%%%%%%%%%%%
disp('方位压缩开始.......');
f_azimuth = fdoc + delta_fdoc + [-Na/2:Na/2-1]/Na*PRF;
d = sqrt(1-(lambda*f_azimuth/2/V).^2);
r0 = distance_range * cos(Theta);
H_azimuth = exp(1i*4*pi/lambda*R0.*d) .* kaiser(Na,2.5)';
signal_final = iftx(signal_rA_src_rcmc .* (H_azimuth.' * ones(1,Nr)));
disp('方位压缩结束');
figure;
% G = 20 * log10(abs(signal_final));                % 换算成分贝（dB）的形式显示
% gm = max(max(G));
% thr06 = 40;                                  % 显示的动态范围40dB,"thr" = threshold
% g_thr = gm - thr06;
% G = (G - g_thr) * (255 / thr06) .* (G > g_thr);
% colormap('gray');
% imagesc(-G);

imagesc(255 - abs(signal_final));
grid on;
title('Final Signal');

%%%%%%%%%%%%%%%%%%%SRC后的距离向PSLR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% H = max(max(abs(signal_rA_src.^2)));
% max_xy = sparse(H == abs(signal_rA_src.^2);
% x = max_xy(1) + [-5:5];
% y = max_xy(2) + [-5:5];

