function f = fourier_ifft_2d(Ft,sx,sy,x,y)
% compute the continuous inverse Fourier transform of function f(x,y):
% Ft(s,xi,eta) = int f(t,x,y) *exp(2*pi*i*xi*x)*exp(2*pi*i*eta*y) dxdy
%
% Algorithm: using the FFT


[Nx,Ny,Nfreq] = size(Ft);

dsx = sx(2)-sx(1);
kx = 0:Nx-1;

dsy = sy(2)-sy(1);
ky = 0:Ny-1;

for n = 1:Nfreq

    for ix = 1:Nx
        for iy = 1:Ny
            Ft(ix,iy,n) = Ft(ix,iy,n).*exp(2*pi*1i*( kx(ix)*dsx*x(1) + ky(iy)*dsy*y(1)));
        end
    end

    f = Nx*Ny*ifftn(Ft); % compute the inverse Discrete Fourier Transform of Ft

    for ix = 1:Nx
        for iy = 1:Ny
            f(ix,iy,n) = f(ix,iy,n)*dsx*dsy.*exp(2*pi*1i*(x(ix)*sx(1) + y(iy)*sy(1)));
        end
    end
end
