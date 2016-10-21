%=================== pulse compression ======================%
%       By Liu Yakun
%                2016.9.20
%============================================================%

clc;clear;close all;

%参数设置  
c = 3e8;   %光速
f0 = 1e9;  %载频1GHz
lambda = c / f0;
Tr = 5e-6;
Br = 100e6;
Kr = Br / Tr;
Fs = 2 * Br;
Ts = 1 / Fs;
N = ceil(Tr / Ts);
Targets = [20,1
           50,1
           80,1
           100,2];  %目标的距离与截面积
R = 200;
t = linspace(-Tr/2,Tr/2,N);
f = linspace(-Fs/2,Fs/2,N);
S_send = exp(1i*pi*Kr*t.^2);

% plot(t,S_send);
S_receive = zeros(1,N);
for i = 1:1:size(Targets,1)
    range = Targets(i,1);
    S_receive = S_receive + Targets(i,2) * exp(-1i*4*pi*range/lambda) .* exp(1i*pi*Kr*(t - 2*range/c).^2);
end

% plot(t,real(S_receive));

h_ref = exp(1i*pi*Kr*t.^2) .* (abs(t) < Tr/2);
H_ref = fft(h_ref);
S_out = ifft(fft(S_receive,N) .* conj(H_ref)) / N;
% plot(t,abs(S_out));
r = Tr * c / 2;
M = ceil(R / r * N);
distance = linspace(0,R,M);

figure;
subplot(311);
plot(t,real(S_send));
subplot(312);
plot(t,real(S_receive));
subplot(313);
plot(distance,abs(S_out(1:M)));
