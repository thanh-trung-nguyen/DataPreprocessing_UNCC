function f = wave_propagation_freq_3d(F,d,freq,x,y,WaveSpeed)
% 
% propagate data in frequency domain: propagate data of the wave equation at z = b
% to z = a. a < b. Here we assume that the wave is outgoing wave propagating in the positive 
% z-direction. 
% INPUT:
% F: data at z = b: F(omega,x,y)
% d: distance of propagation in METERS: d = b-a > 0
% freq: frequency in Hz (a vector) or GHz, depending on unit used. If it is
% in GHz, the wave speed must be in m/ns. 
% x, y: coordinates of the data points in the x and y directions (in METERS)
% WaveSpeed: wave speed (m/s or m/ns).
% OUTPUT: f: the propagated dat at x = a. 


[Nfreq,Nx,Ny] = size(F);

if Nx ~=length(x)
    error('The second dimension of F must be equal to the length of x');
end
if Ny~=length(y)
    error('The third dimension of F must be equal to the length of y');
end


% compute the Fourier transform of the function F at z = b
FF = 0*F; 
for n = 1:Nfreq
    [Fourier,sx,sy] = fourier_fft_2d(reshape(F(n,:,:),Nx,Ny),x,y);
    FF(n,:,:) = reshape(Fourier,1,Nx,Ny);
end

FFF = 0*FF;

kSquared = (freq./WaveSpeed).^2; % square of (wavenumber/2pi)

for idy = 1:Ny
    for idx = 1:Nx
        v = kSquared  - sx(idx)^2 - sy(idy)^2;
        idf = find(v > 0);        
        if ~isempty(idf)
            FFF(idf,idx,idy) = FF(idf,idx,idy).*exp(2*pi*1i*d*sqrt(v(idf)));   
        end
    end
end

f = 0*FFF;
for n = 1:Nfreq
    Fourier = fourier_ifft_2d(reshape(FFF(n,:,:),Nx,Ny),sx,sy,x,y);
    f(n,:,:) = reshape(Fourier,1,Nx,Ny);
end




