function f = fourier_ifft_3d(Ft,st,sx,sy,t,x,y)
% compute the continuous inverse Fourier transform of function f(t,x,y):
% Ft(s,xi,eta) = int f(t,x,y) exp(2*pi*i*s*t)*exp(2*pi*i*xi*x)*exp(2*pi*i*eta*y) dtdxdy
%
% Algorithm: using the FFT


[Nt,Nx,Ny] = size(Ft);
dst = st(2)-st(1);
kt = 0:Nt-1;

dsx = sx(2)-sx(1);
kx = 0:Nx-1;

dsy = sy(2)-sy(1);
ky = 0:Ny-1;

if size(t,1) == 1
    t = t';
end

for ix = 1:Nx
    for iy = 1:Ny
        Ft(:,ix,iy) = Ft(:,ix,iy).*exp(2*pi*1i*(kt.'*dst*t(1) + kx(ix)*dsx*x(1) + ky(iy)*dsy*y(1)));
    end
end

f = Nt*Nx*Ny*ifftn(Ft); % compute the inverse Discrete Fourier Transform of Ft

for ix = 1:Nx
    for iy = 1:Ny
        f(:,ix,iy) = f(:,ix,iy)*dst*dsx*dsy.*exp(2*pi*1i*(t*st(1) + x(ix)*sx(1) + y(iy)*sy(1)));
    end
end
