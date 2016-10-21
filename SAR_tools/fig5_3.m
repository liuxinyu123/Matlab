% use this program to reproduce Fig. 5.2 of text
clear all
close all
nscat = 2; %two point scatterers
taup = 10e-6; % 100 microsecond uncompressed pulse
b = 50.0e6; % 50 MHz bandwdith
rrec = 50 ; % 50 meter processing window
scat_range = [15 25] ; % scattterers are 15 and 25 meters into window
scat_rcs = [1 2]; % RCS 1 m^2 and 2m^2
winid = 0; %no window used
[y] = matched_filter(nscat,taup,b,rrec,scat_range,scat_rcs,winid);