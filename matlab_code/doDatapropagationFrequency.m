% data propagation of the frequency domain data: 
close all; clear;


% -----------the following parameters can be changed.
% datafile = 'bub-bam1-lateral.csv'; % you can change this file name to try another data set
% datafile = 'waterbottle-in-air.csv'; % you can change this file name to try another data set
datafile = 'bamboo1-in-air.csv';

emptyspacefile = 

PropDistance = 0.4; % in METERS, you may change this. 
FreqList = [2.09e9, 2.45e9, 2.91e9, 4.55e9, 7e9]; % choose indices of frequencies to be displayed for the above distance of propagation

Nx = 50; % number of grid points in x
Ny = 51; % number of grid point in y
Xmin = -0.5;
Xmax = 0.5; 
Ymin = -0.5; 
Ymax = 0.5; 
Nfreq = 100;

ListOfDistance = 0.1:0.1:1; % test with different distances. 
FreqTest = 7e9; % frequency at which the data is propagated at different distances

% -----------


LightSpeed = 299792458; % m/s.

if ~exist(datafile,'file');
    error('Data file is not found. You must move to the correct folder');
end


% load the data into Matlab:  
[freq,data3d,data2d] = load_data_freq(datafile,Nfreq,Nx,Ny);

X = linspace(Xmin,Xmax,Nx); Y = linspace(Ymin,Ymax,Ny);

dataProp = wave_propagation_freq_3d(data3d,PropDistance,freq',X,Y,LightSpeed);
dataProp = reshape(dataProp,Nfreq,Nx*Ny);


% ================================================
% display the result: 

% figure; 
% for n = 1:Nfreq; 
%     imagesc(reshape(real(dataProp(n,:,:)),Nx,Ny)'); 
%     title(freq(n)); colorbar; 
%     pause(0.1); 
% end

% -------------save some figures: 
figure; set(gca,'fontsize',15); 
figname = datafile(1:end-4);


% the original data: 
for k = 1:length(FreqList);
    x = abs(freq - FreqList(k));
    n = find(x ==min(x));
    imagesc(real(reshape(data3d(n,:,:),Nx,Ny)')); 
    title(sprintf('%s%5.2f%s','Real part. Frequency = ',freq(n)/10^9,' GHz'));
    colorbar; 
    xlabel('x'); ylabel('y'); 
    print('-depsc2',sprintf('%s%s%4.2f%s',figname,'_origdat_freq',freq(n)/10^9,'.eps'));
end


% propagated data at a fixed distance of propagation:
for k = 1:length(FreqList);
    x = abs(freq - FreqList(k));
    n = find(x ==min(x));
    imagesc(real(reshape(dataProp(n,:,:),Nx,Ny)')); 
    title(sprintf('%s%5.2f%s%4.2f%s','Real part. Frequency = ',freq(n)/10^9,' GHz. Propagation distance = ',PropDistance,' (m)'));
    colorbar; 
    xlabel('x'); ylabel('y'); 
    print('-depsc2',sprintf('%s%s%4.2f%s%4.2f%s',figname,'_datprop_dist',PropDistance','_freq',freq(n)/10^9,'.eps'));
end



% ========= test with different propagation distances:
x = abs(freq - FreqTest); % 7GHz
FreqIdx = find(x==min(x));

for k = 1:length(ListOfDistance)
    PropDistance = ListOfDistance(k);
    dataProp = wave_propagation_freq_3d(data3d(FreqIdx,:,:),PropDistance,freq(FreqIdx),X,Y,LightSpeed);

    imagesc(real(reshape(dataProp,Nx,Ny)')); 
    title(sprintf('%s%5.2f%s%4.2f%s','Real part. Frequency = ',freq(FreqIdx)/10^9,' GHz. Propagation distance = ',PropDistance,' (m)'));
    colorbar; 
    xlabel('x'); ylabel('y'); 
    print('-depsc2',sprintf('%s%s%4.2f%s%4.2f%s',figname,'_datprop_dist',PropDistance','_freq',freq(FreqIdx)/10^9,'.eps'));
    disp('Press any key to continue');
    pause;
    
end
