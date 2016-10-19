function y = fty(x)
y = fftshift(fft(fftshift(x')))';