function [Ft, s] = FourierTransformFFT(F,t,s_min,Type,ds)
% FourierTransformFFT
% calculate the 1d continuous Fourier transform
%     Ft(s) = int_{\-infty}^{\infty} F(t) exp(-2*pi*t*s) dt
% or the inverse Fourier transform: 
%     Ft(s) = int_{\-infty}^{\infty} F(t) exp(2*pi*t*s) dt
% 
% Algorithm: use FFT or IFFT
% 
% Input: 
% 1. F: a vector, 2D or 3D matrix of the original function F(t). If F is a matrix, the Fourier transform is calculated for the first index
% 2. t: a vector of variable t with a constant step size. Note that t may mean time or frequency.
% 3. s_min: the minimum value of the variable after Fourier transform, other values of s are determined by FFT. More precisely, 
%     s = s_min + 0:1:(Nt-1)*ds, with ds = 1/Nt/dt, where Nt is the number of data points w.r.t. t, dt: step size in t. 
% 4. Type: if Type = -1, Ft is the inverse Fourier transform (the second formula). Otherwise, the first formula is used 
% 5. ds (optional): the desired step size in output variable s. Note that the FFT algorithm calcualtes s automatically. 
% However, the step size in s can be reduced by padding zeros at the end of the input data F. If this parameter is larger
% than 1/N/dt, the latter is taken. 
% 
% OUTPUT: 
% Ft: vector or matrix of Fourier transform of the same size as F. 
% s: vector of s values at which the FT is calculated. 
% 
% CALL THE FUNCTION: 
%   [Ft, s] = FourierTransformFFT(F,t,s_min): use the first formula
%   [Ft, s] = FourierTransformFFT(F,t,s_min,Type)
%   [Ft, s] = FourierTransformFFT(F,t,s_min,Type,ds)
%
% NOTE: if s_min = 0, only the data of the first half of the frequencies
% are usable because Ft(Nt-k+2) = conj(Ft(k)), k >=2 from the FFT algorithm
%  This function can replace the routines "fourier_fft" and "fourier_ifft".
%  
%
% ======= Nguyen Trung Thanh, 2016. 

[N1,N2,N3] = size(F); % F could be a 3D matrix.
Nt = length(t); 

% check for consistency of the input parameters:
if Nt ~= N1
    error('Error in FourierTransformFFT: The length of t must be the same as the number of rows of F');
end

dt = t(2)-t(1); % time step

% default Type: forward Fourier transform: 
if nargin < 4
    Type = 1;
end

% If ds is smaller than the default stepsize, add zeros to the end of the input data.
if (nargin == 5) && (ds < 1/(Nt*dt))
  Nt = ceil(1/(ds*dt));
  F = [F; zeros(Nt-N1,N2,N3)];
end

n = (0:1:Nt-1)'; % n is a column vector 

% compute the frequencies:
ds = 1/(Nt*dt); 
s = n*ds + s_min; % vector of frequencies

% compute the Fourier transform:
Ft = 0*F; 
if Type == -1 % inverse Fourier transform: 
    for n3 = 1:N3
      F(:,:,n3) = F(:,:,n3).*(exp(2*pi*1i*s_min*n*dt)*ones(1,N2));
      Ft(:,:,n3) = ifft(F(:,:,n3)); 
      Ft(:,:,n3) = Ft(:,:,n3).*dt.*Nt.*(exp(2*pi*1i*t(1)*s)*ones(1,N2));
    end    
else    
    for n3 = 1:N3
      F(:,:,n3) = F(:,:,n3).*(exp(-2*pi*1i*s_min*n*dt)*ones(1,N2));
      Ft(:,:,n3) = fft(F(:,:,n3)); 
      Ft(:,:,n3) = Ft(:,:,n3)*dt.*(exp(-2*pi*1i*t(1)*s)*ones(1,N2));
    end
end
