% 2016/10/12

clc;
clear;
close all;
% 雷达平台参数
C = 3e8;
Fc = 5.3e9;
Vr = 150;
theta_rc = 0 / 180 * pi; %波束斜视角
H = 5000;
R0 = 2e4;
Y0 = 1e4;%场景中心Y轴
La = 4;


% 快时间参数
Tr = 2.5e-6;
Kr = 20e12;
alpha_Fsr = 1.2;%  距离过采样率

% 慢时间参数
alpha_Fsa = 1.25;
delta_x = 400;
delta_y = 150;

Targets = [ 0  -100   1
            200 0     2
            0   100   2
           -200 0  1];



