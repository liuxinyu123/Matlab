function S = echo_creation(C,H,Yc,lambda,Lsar,Kr,Tr,Tau,Ra,Targets)

%C 光速
%H  雷达飞行高度
%Yc  场景中心Y轴坐标
%lambda  发射信号波长
%Lsar 合成孔径长度
%Kr  发射信号调频率
%Tr 脉冲持续时间
%Tau 距离向时间向量
%Ra 方位向距离向量
%Targets 目标参数 矩阵 格式为 行： x ，y，rcs  ；列：每个目标

N = size(Targets,1);%目标数量
Nr = size(Tau,2); %距离向采样点数
Na = size(Ra,2); %方位向采样点数

S = zeros(Na,Nr); %初始化回波数据
for i = 1:N
    delta_x = Ra - Targets(i,1); %平台距目标的x轴距离差
    delta_y = Targets(i,2);%平台距目标的y轴距离差
    range = sqrt(delta_x.^2 + H^2 + delta_y^2); % 每个方位采样点 雷达到目标的斜距
    rcs = Targets(i,3); %点目标的后向散射系数
    tau = ones(Na,1) * Tau - (2*range / C)' * ones(1,Nr); % 时间差矩阵
    phase = pi * Kr * tau.^2 - 4 * pi / lambda * (range' * ones(1,Nr));%接收信号相位
%     S = S + rcs * exp(1i * phase) .* ((abs(delta_x) < Lsar / 2)' * ones(1,Nr)) .* (tau > 0 & tau < Tr);
    S = S + rcs * exp(1i * phase) .* ((abs(delta_x) < Lsar / 2)' * ones(1,Nr)) .* (abs(tau) < Tr/2);
end

