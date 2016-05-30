function [dataFreq,freq] = dataTime2Freq(dataTime,Time)



[Nt,Nx] = size(dataTime);

dataFreq = 0*dataTime; 

for n = 1:Nx
    [Ft,freq] = fourier_fft(dataTime(:,n),Time);
    dataFreq(:,n) = reshape(Ft,Nt,1);
end

Nfreq = floor((Nt-1)/2);

dataFreq = dataFreq(Nt-Nfreq+1:Nt,:);
freq = freq(Nt-Nfreq+1:Nt);