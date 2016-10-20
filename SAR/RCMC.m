function y = RCMC(signal_rD,C,lambda,V,R0,Fs,PRF,mode)

%signal_rD  RCMC输入信号 为距离多普勒表示
%C 光速
%lambda 脉冲波长
%V 雷达飞行速度
%R0 雷达航线到场景中心的最短斜距
%Fs 距离向采样率
%PRF 脉冲重复发射频率
%DY 距离向分辨率
%mode 插值模式  1 为最近邻域插值   2 为sinc插值
y = signal_rD;

Na = size(y,1);
Nr = size(y,2);

if mode == 1
    win = waitbar(0,'最近邻域插值');
    for i = 1:Na
        for j = 1:Nr
            delta_R = (1/8)*(lambda/V)^2*(R0+(j-Nr/2)*C/2/Fs)*((j-Nr/2)/Nr*PRF)^2;%距离徙动量
            RCM = delta_R * 2 * Fs / C;%徙动了多少个距离单元
            delta_RCM = RCM - floor(RCM);%小数部分
            if round(RCM + j) > Nr
                y(i,j) = y(i,Nr/2);
            else
                if delta_RCM < 0.5
                    y(i,j) = y(i,j+floor(RCM));
                else
                    y(i,j) = y(i,j+ceil(RCM));
                end
            end
        end
        waitbar(i/Na);
    end
    close(win);
end

if mode == 2
    win = waitbar(0,'sinc插值');
    N_core = 4;
    y = zeros(Na,Nr);
    for i = 1:Na
        for j = 1:Nr
            delta_R = 1/8*(lambda/V)^2*(R0+(j-Nr/2)*C/2/Fs)*((j-Nr/2)/Nr*PRF)^2;%距离徙动量
            RCM = delta_R / DY;%徙动了多少个距离单元
            delta_RCM = RCM - floor(RCM);%小数部分
            for k = -N_core/2 : N_core/2-1
                if j+k+RCM > Nr
                    y(i,j) = y(i,j) + signal_rD(i,j)*sinc(k+RCM);
                else
                    y(i,j) = y(i,j) + signal_rD(i,j+floor(RCM)+k)*sinc(k+delta_RCM);
                end
            end
        end
        waitbar(i/Na);
    end
    close(win);
end

    