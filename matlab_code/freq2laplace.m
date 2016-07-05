function u = freq2laplace(FreqData,Freq,s)
% function FREQ2LAPLACE 
% convert from Frequency domain data to Laplace transformed data. 
% The relationship between the frequency data and the Laplace transformed
% data is: 
% u(x,s) = 1/pi*real(int_0^infty uf(x,k)/(s + ik)dk)
% 
% Input: 
% FreqData: a 2d or 3d matrix. The first index (column) is the frequency data
% Freq: a vector of the frequencies of the data
% s: a vector of pseudo frequency in the Laplace transform. 
% 
% Output: 
% u: a matrix (2d or 3d, depending on the input matrix) with the first index 
% being the Laplace transform at each location. 
% 


Ns = length(s); % number of pseudo frequencies
Nf = length(Freq); % number of real frequencies

df = Freq(2) - Freq(1); % frequency step
Freq = Freq(:); % convert to column vector if needed. 

[Nf2,Nx,Ny] = size(FreqData);

if Nf2 ~=Nf
    error('Number of frequencies must be equal to the first size of the input data matrix');
end

u = zeros(Ns,Nx,Ny); 

for ids = 1:Ns
    for iy = 1:Ny
        for ix = 1:Nx
            u(ids,ix,iy) = 2*real(sum(FreqData(:,ix,iy)./(s(ids) - 2*pi*1i*Freq))*df);
        end
    end
end


