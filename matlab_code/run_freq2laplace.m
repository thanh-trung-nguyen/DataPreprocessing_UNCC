Nx = 50; Ny = 51; Nfreq = 100; 
Xmin = -0.5; Xmax = 0.5; 
Ymin = -0.5; Ymax = 0.5; 
X = linspace(Xmin,Xmax,Nx); 
Y = linspace(Ymin,Ymax,Ny);

datafolder = '/home/thanhnguyen/Programming_Git/UNCC_data/FreqData_2016/non-blind/Objects-in-air/Data-files/';
filename = 'waterbottle-in-air.csv';
[freq,data3d,data2d] = load_data_freq([datafolder,filename],100,50,51);
PropDist = 0.9; % (m)
LightSpeed = 299792458; % m/s.


dataProp = wave_propagation_freq_3d(data3d,PropDist,freq',X,Y,LightSpeed);
dataProp = reshape(dataProp,Nfreq,Nx*Ny);
% convert to Laplace transform domain: 
s = 8; 
us = freq2laplace(dataProp,freq/10^9,s); 
us(us<0) = 0; 
imagesc(reshape(us,Nx,Ny)'); colorbar;  title(sprintf('%s%d','s = ',s)); 
print('-dpng',[filename(1:end-4),'_laplace.png']);
