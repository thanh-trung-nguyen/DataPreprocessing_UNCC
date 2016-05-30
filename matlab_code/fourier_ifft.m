function f = fourier_ifft(Ft,s,t)
% compute the continuous inverse Fourier transform of function f:
% f(t) = int_\-infty^\infty Ft(s) exp(2*pi*i*s*t)ds
% for t_n = t
% Algorithm: using the FFT
% output: a ROW vector f 


if size(Ft,2) == 1
    Ft = Ft.';
end

if size(t,2) == 1
    t = t.';
end
if size(s,2) == 1
    s = s.';
end


Nt = length(Ft);
if Nt ~= length(t)
    error('Time variable must have the same length as the frequency variable');
end

ds = s(2)-s(1);
if (t(1)~=0)
    k = 0:Nt-1;
    Ft = Ft.*exp(2*pi*1i*k*ds*t(1));
end
f = Nt*ifft(Ft); % compute the inverse Discrete Fourier Transform of Ft
f = f*ds.*exp(2*pi*1i*t*s(1));







