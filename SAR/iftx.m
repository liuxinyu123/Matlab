function y = iftx(x)
y = fftshift(ifft(fftshift(x)));