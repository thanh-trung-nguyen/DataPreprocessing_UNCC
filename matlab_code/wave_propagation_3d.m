function f = wave_propagation_3d(F,d,t,x,y,WaveSpeed)
% WAVE_PROPAGATION_3d
% propagates data from one plane to another plane. Note that the propagation is 
% in the opposite direction as the direction of wave propagation.
% F: a 3D matrix. The first index is for time, the second index is for x, the last index is the y direction.
% d: distance of propagation
% t: vector of time instant
% x, y: vectors of grid points in the x and y directions
% WaveSpeed (optional): Speed of wave in the medium. If this parameter is not provided, it is assumed to be 1 (normalized).  
% How to call: 
% f = wave_propagation_3d(F,d,t,x,y). 
% f = wave_propagation_3d(F,d,t,x,y,WaveSpeed)
% 
% ==========================



if (nargin < 6)
    WaveSpeed = 1;
end

[Nt,Nx,Ny] = size(F);

% compute the Fourier transform of the function F at z = d
[FF,st,sx,sy] = fourier_fft_3d(F,t,x,y);

if size(st,1)==1
    st = st.';
end
stSquared = (st/WaveSpeed).^2;


% idx = 1:length(st);
% compute Fjkl = F^jkl*exp(2pi*i*a*sqrt()):

FFF = 0*FF;
for xx = 1:Nx
    for yy = 1:Ny
        v = stSquared  - sx(xx)^2 - sy(yy)^2;
        idx = find(v > 0);        
        if ~isempty(idx)
            FFF(idx,xx,yy) = FF(idx,xx,yy).*exp(2*pi*1i*d*sqrt(v(idx)));  
        end
    end
end

f = real(fourier_ifft_3d(FFF,st,sx,sy,t,x,y));










