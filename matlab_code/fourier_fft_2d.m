function [f,fre1,fre2] = fourier_fft_2d(F,x,y)
% compute the continuous Fourier transform of function f(x,y):
% Ft(xi,eta) = int f(x,y) *exp(-2*pi*i*xi*x)*exp(-2*pi*i*eta*y) dxdy
% INPUT:
% F: the first index is for x, the second index is for y. Moreover, F may depend on a third index which
% is used for frequencies.
%
% Algorithm: using the FFT


[Nx,Ny,Nfreq] = size(F);
dx = x(2)-x(1);
dy = y(2)-y(1);


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
f = zeros(Nx,Ny,Nfreq);
for n = 1:Nfreq
  for ix = 1:Nx
      for iy = 1:Ny
          F(ix,iy,n) = F(ix,iy,n).*exp(-2*pi*1i*(sx_min*nx(ix)*dx + sy_min*ny(iy)*dy));
      end
  end

  f(:,:,n) = fftn(F(:,:,n));

  for ix = 1:Nx
      for iy = 1:Ny
          f(ix,iy,n) = f(ix,iy,n)*dx*dy.*exp(-2*pi*1i*(sx(ix)*x(1) + sy(iy)*y(1)));
      end
  end
end

if nargin > 1
    fre1 = sx;
    fre2 = sy;
end







