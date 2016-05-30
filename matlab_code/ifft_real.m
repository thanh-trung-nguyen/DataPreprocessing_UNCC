function [F,t] = ifft_real(Fs,s,t_min)
% FFT_REAL
% computes the inverse Fourier transform of the form: 
% F(t) = int_{\-infty}^{\infty} Fs(s) exp(2*pi*t*s) ds
% assuming that F(t) is a real function. In other words, Fs is the complex function which is the Fourier transform of the
% real function F(t). 
% 
% Input: 
% 1. F: a vector, 2D or 3D matrix of the original function F(t). If F is a matrix, the Fourier transform is calculated for the first index
% 2. s: a vector of frequencies with a constant step size. Note that s is positive.
% 3. t_min: the minimum time. Default is zero. 
% 
% OUTPUT: 
% F: vector or matrix of inverse Fourier transform. 
% t: vector of t values at which the FT is calculated. 
% 
% CALL THE FUNCTION: 
%   [Ft, s] = ifft_real(Fs,s): use the first formula
%   [Ft, s] = ifft_real(Fs,s,t_min)
%
% Note: Fs is only the measured data for positive frequencies. Since F is
% assumed to be real, F(-s) = conj(F(s)). Therefore, we must expand F to
% negative frequencies as well before calculating the inverse Fourier
% transform. Moreover, if the data is given for s in (s_min, s_max), where
% s_min > 0, then we must assume that F(s) = 0 for s < s_min, and expand
% the data to frequencies s, 0 < s < s_min by adding zeros. 
%
% ======= Nguyen Trung Thanh, 2016. 


% expand the frequency domain data: 
ds = s(2) - s(1); 
if s(1) > ds
    s = [zeros(
