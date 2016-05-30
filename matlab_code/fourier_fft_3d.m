function [f,fre1,fre2,fre3] = fourier_fft_3d(F,t,x,y)
% compute the continuous Fourier transform of function f(t,x,y):
% Ft(s,xi,eta) = int f(t,x,y) exp(-2*pi*i*s*t)*exp(-2*pi*i*xi*x)*exp(-2*pi*i*eta*y) dtdxdy
%
% Algorithm: using the FFT


[Nt,Nx,Ny] = size(F);
dt = t(2)-t(1);
dx = x(2)-x(1);
dy = y(2)-y(1);


% compute the frequencies in t:
nt = 0:1:Nt-1;
dst = 1/Nt/dt; 
st_min = -(floor((Nt-1)/2))*dst;
st = nt*dst + st_min;

% compute the frequencies in x:
nx = 0:1:Nx-1;
dsx = 1/Nx/dx; 
sx_min = -(floor((Nx-1)/2))*dsx;
sx = nx*dsx + sx_min;

% compute the frequencies in y:
ny = 0:1:Ny-1;
dsy = 1/Ny/dy; 
sy_min = -(floor((Ny-1)/2))*dsy;
sy = ny*dsy + sy_min;


% compute the Fourier transform:
for ix = 1:Nx
    for iy = 1:Ny
        F(:,ix,iy) = F(:,ix,iy).*exp(-2*pi*1i*(st_min*nt.'*dt + sx_min*nx(ix)*dx + sy_min*ny(iy)*dy));
    end
end

f = fftn(F);

for ix = 1:Nx
    for iy = 1:Ny
        f(:,ix,iy) = f(:,ix,iy)*dt*dx*dy.*exp(-2*pi*1i*(st.'*t(1) + sx(ix)*x(1) + sy(iy)*y(1)));
    end
end


if nargin > 1
    fre1 = st;
    fre2 = sx;
    fre3 = sy;
end







