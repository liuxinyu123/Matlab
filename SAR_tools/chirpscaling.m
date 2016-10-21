% ========================================================================
% Chirp Scaling Algorithm
% ========================================================================
% Implementation of the basic algorithm for Extended Chirp Scaling after
% the paper:
% "Extended Chirp Scaling Algorithm for Air- and Spaceborne SAR Data
% Processing in Stripmap and ScanSAR Imaging Modes". A. Moreira,
% J. Mittermayer, R. Scheiber. IEEE Trans Geoscience and Remote Sensing,
% 34(5): 1123-1136, 1996
% ========================================================================

clear all; 
close all;
doPlot = 0; % 1 to draw plots, 0 to skip all plots
doMesh = 0; % 1 to draw meshs, 0 to skip all meshs

% ========================================================================
% set parameters, read parameters
% ========================================================================
% example raw data point target from file 'pointtarget.raw'
c = 299702547; % in air, vacuum: 299792458;
f_c = 10000000000; % carrier frequency
f_s = 100000000; % data sampling rate
echoes = 348; % number of echoes in data file
samples = 152; % number of samples per echo
tau = 2.83529480918429e-006:1/f_s:2.83529480918429e-006+samples/f_s-1/f_s; % fast time
f_dc = 0; % Doppler centroid
v = 200; % SAR platform velocity
PRF = 1000; % pulse repetition frequency
t_p = 5e-7; % chirp pulse duration
B = 100000000; % chirp bandwidth

% envisat swath IS7 IM Mode:
% f_c = 5331004416;
% f_s = 19207680;
% echoes = 2048;
% samples = 2048;
% tau = 0.00691589:1/f_s:0.00691589+samples/f_s-1/f_s;
% f_dc = -97.341634;
% v = 7065.51587;
% PRF = 2067.120103315;
% t_p = 2.160594095e-5;
% B = 16000000;

r_ref = (tau(1)+samples/2/f_s)/2*c;%--------------参考距离
alpha = 1;
f_a = -PRF/2+f_dc:PRF/echoes:f_dc+PRF/2-PRF/echoes;
f_r = -f_s/2:f_s/samples:f_s/2-f_s/samples;
lambda = c/f_c;

% read the raw data
data=readMatrix('pointtarget.raw',1,samples,1,echoes);
data0=data;
G=abs(data);
xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);
figure;
colormap(gray(256))
image(256-cg*(G-ng));
axis('image');axis('xy')
xlabel('range')
ylabel('azimuth')
title('SAR raw data')

% ========================================================================
% azimuth fft
% ========================================================================
data = ftx(data);

if (doPlot==1)
    G=angle(data);
    xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);
    figure;
    colormap(jet(256))
    image(256-cg*(G-ng));
    axis('image');axis('xy')
    xlabel('range')
    ylabel('Doppler frequency')
    title('range signal/Doppler domain')
end;
%%

% ========================================================================
% chirp scaling, range scaling: H1
% ========================================================================
beta = (1 - (f_a*lambda/2/v).^2).^0.5;
a = 1./beta - 1;%----（4-29式）
R = r_ref./beta;%----（4-28式）
a_scl = a + (1-alpha).*(1+a)./alpha;
k_r = -B./t_p;
% WARNING: in the work of Moreira et al, 1996 k_r is defined negative
% because they assume a "down-chirp", whereas we assume an "up-chirp" and
% must explicitely set k_r as negative for the following transfer functions
% to be equal to the theoretical formulations.

k_inv = 1./k_r - (2.*lambda.*r_ref.*(beta.^2-1))./(c^2.*beta.^3);%----（4-31式）
k = 1./k_inv;
%k = (k_r*c.^2.*beta.^3) ./ (c.^2.*beta.^3 - ...
%k_r.*2.*lambda.*r_ref.*(beta.^2-1)); % same thing as above

x = k.*a_scl;
X = x(:)*ones(1,samples);
Tau = ones(echoes,1)*tau;
y = 2.*R./c;
Y = y(:)*ones(1,samples);
Z = (Tau-Y).^2;
H1 = exp(-j*pi*X.*Z);

data = data .* H1;%完成第一次相位相乘
%%

if (doPlot==1)
    G=angle(data);
    xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);
    figure;
    colormap(jet(256))
    image(256-cg*(G-ng));
    axis('image');axis('xy')
    xlabel('range')
    ylabel('Doppler frequency')
    title('chirp scaling and range scaling')
end;

if (doMesh==1)
    figure;
    colormap(jet(256))
    mesh(abs(iftx(data)));
    xlabel('range')
    ylabel('azimuth')
    title('H1')
end;

% ========================================================================
% range fft
% ========================================================================
data = fty(data);

if (doPlot==1)
    G=angle(data);
    xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);
    figure;
    colormap(jet(256))
    image(256-cg*(G-ng));
    axis('image');axis('xy')
    xlabel('signal frequency')
    ylabel('Doppler frequency')
    title('2D spectral domain')
end;

% ========================================================================
% bulk rcmc, range compression: H2
% ========================================================================
x = 1./(k.*(1+a_scl));
X = x(:)*ones(1,samples);
y = f_r.^2;
Y = ones(echoes,1)*y;
z = f_r;
Z = ones(echoes,1)*z;
A = a(:)*ones(1,samples);
Z = Z .* A;
H2 = exp(-j*pi*X.*Y) .* exp(j*4*pi*r_ref/c.*Z);
data = data .* H2;

if (doPlot==1)
    G=angle(data);
    xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);
    figure;
    colormap(jet(256))
    image(256-cg*(G-ng));
    axis('image');axis('xy')
    xlabel('signal frequency')
    ylabel('Doppler frequency')
    title('bulk rcmc and range compression')
end;

if (doMesh==1)
    figure;
    colormap(jet(256))
    mesh(abs(iftx(ifty(data))));
    xlabel('range')
    ylabel('azimuth')
    title('H2')
end;

% ========================================================================
% range ifft
% ========================================================================
data = ifty(data);

if (doPlot==1)
    G=angle(data);
    xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);
    figure;
    colormap(jet(256))
    image(256-cg*(G-ng));
    axis('image');axis('xy')
    xlabel('range')
    ylabel('Doppler frequency')
    title('range Doppler Domain after range compression')
end;
%%

% ========================================================================
% angle correction: H3
% ========================================================================
r_0 = tau/2*c;
x = k.*a_scl.*(1+a).^2 ./ (c^2.*(1+a_scl));
X = x(:)*ones(1,samples);
z = (r_0-r_ref).^2;
Z = ones(echoes,1)*z;
dphi = 4*pi*X.*Z;
H3 = exp(j*dphi);

data = data .* H3;

if (doPlot==1)
    G=angle(data);
    xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);
    figure;
    colormap(jet(256))
    image(256-cg*(G-ng));
    axis('image');axis('xy')
    xlabel('range')
    ylabel('Doppler frequency')
    title('range Doppler domain after angle correction')
end;

if (doMesh==1)
    figure;
    colormap(jet(256))
    mesh(abs(iftx(data)));
    xlabel('range')
    ylabel('azimuth')
    title('H3')
end;

% ========================================================================
% azimuth compression: H4
% ========================================================================
r_0_scl = r_ref + (r_0-r_ref)/alpha;
X = ones(echoes,1)*r_0_scl;
Z = (beta(:)-1)*ones(1,samples);
H4 = exp(j*4*pi/lambda*X.*Z);

data = data .* H4;

if (doPlot==1)
    G=angle(data);
    xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);
    figure;
    colormap(jet(256))
    image(256-cg*(G-ng));
    axis('image');axis('xy')
    xlabel('range')
    ylabel('Doppler frequency')
    title('azimuth compression')
end;

if (doMesh==1)
    figure;
    colormap(jet(256))
    mesh(abs(iftx(data)));
    xlabel('range')
    ylabel('azimuth')
    title('H4')
end;
% ========================================================================
% azimuth ifft
% ========================================================================
data = iftx(data);

% ========================================================================
% write result
% ========================================================================

G=abs(data);
xg=max(max(G)); ng=min(min(G)); cg=255/(xg-ng);
figure;
% colormap(gray(256))
% image(256-cg*(G-ng));
mesh(G);
axis('image');axis('xy')
xlabel('range')
ylabel('azimuth')
title('focused data')

% = end ==================================================================
