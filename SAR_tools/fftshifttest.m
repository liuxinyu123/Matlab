clear;

clc;
close all;
t=0:0.001:2;

n=2001;

Fs=1000;

Fc=200;

x=cos(2*pi*Fc*t);

y1=fft(x);

y2=fftshift(y1);

y3=fftshift(x);

y4=fft(y3);

f=(0:2000)*Fs/n-Fs/2;

hold on;

plot(f,abs(y1),'r'); 

plot(f,abs(y2),'b');

figure;
hold on;
plot(f,abs(y3),'g');

plot(f,abs(y4),'y');