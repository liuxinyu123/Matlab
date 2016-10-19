function y = pulse_compression(signal_org,signal_ref,Na)

% signal_org  接收到的信号
%signal_ref 参考的信号
%方位向采样点数

signal_Ra = fty(signal_org);
signal_REF = fty(signal_ref);

y = ifty(signal_Ra .* conj(ones(Na,1) * signal_REF));
