function dat_time = dataFreq2Time(dataFreq,Freq,Time)
% convert from frequency domain data to time domain data





Npoints = size(dataFreq,1); % number of measurement points.
Ntime = length(Time);

dat_time = zeros(Npoints,Ntime);

for n = 1:Npoints
    dat_time(n,:) = fourier_ifft(dataFreq(n,:),Freq,Time);
end


