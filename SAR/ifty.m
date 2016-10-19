function y = ifty(x)
y = fftshift(fft(fftshift(x')))';