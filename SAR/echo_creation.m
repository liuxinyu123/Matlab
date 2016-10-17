function S = echo_creation(C,H,Yc,lambda,Lsar,Kr,Tr,Tau,Ra,Targets)
%Lsar 合成孔径长度
%Tr 脉冲持续时间
%Tau 距离向时间向量
%Ra 方位向距离向量
%Targets 目标参数 矩阵 格式为 行： x ，y，rcs  ；列：每个目标
N = size(Targets,1);%目标数量
Nr = size(Tau,2);
Na = size(Ra,2);

S = zeros(Na,Nr);
for i = 1:N
    delta_x = Ra - Targets(i,1); %平台距目标的x轴距离差
    delta_y = Yc;%y轴距离差
    delta_z = H;%z轴距离差
    
    range = sqrt(delta_x.^2 + delta_y^2 + delta_z^2);
    rcs = Targets(i,3);
    tau = ones(Na,1) * Tau - (2*range ./ C)' * ones(1,Nr); % 时间差矩阵
    phase = pi * (-4 / lambda * (range' * ones(1,Nr)) + Kr * tau.^2);%接收信号相位
    S = S + rcs * exp(1i * phase) .* ((abs(delta_x) < Lsar / 2)' * ones(1,Nr)) .* (tau > 0 & tau < Tr);
end

