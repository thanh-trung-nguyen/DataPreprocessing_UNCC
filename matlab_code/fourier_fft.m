function [f,freq] = fourier_fft(F,t)
% compute the continuous Fourier transform of function f:
%Ft(s) = int_\-infty^\infty f(t) exp(-2*pi*i*s*t)dt or its inverse
%
% Algorithm: using the FFT
% output: a ROW vector f and a ROW vector of frequencies


% convert to row vectors if needed:
if size(F,2) == 1
    F = F.';
end
if size(t,2) == 1
    t = t.';
end
Nt = length(t);

if Nt ~= length(F)
    error('Vectors of variable and function value must be of the same length');
end

dt = t(2)-t(1);

n = 0:1:Nt-1;

% compute the frequencies:
ds = 1/Nt/dt; 
s_min = -(floor((Nt-1)/2))*ds; % zero frequency is in the centre of the frequency vector
% s_min = 0;
st = n*ds + s_min;


% compute the Fourier transform:
F = F.*exp(-2*pi*1i*s_min*n*dt);
f = fft(F); 
f = f*dt.*exp(-2*pi*1i*t(1)*st);


if nargin > 1
    freq = st;
end







