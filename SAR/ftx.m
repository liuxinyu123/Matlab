function y = ftx(x)
y = fftshift(fft(fftshift(x)));
