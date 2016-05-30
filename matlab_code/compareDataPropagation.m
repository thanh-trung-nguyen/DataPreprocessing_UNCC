% compare data propagation in time domain and frequency domain:
% Test for a time-dependent data set: pre-processed-air-2. Before you
% can run this script, you must go to this folder "pre-processed-air-2".
% You can test with other data sets too, just change the name of the
% datafile. 



datafile = 'object_8/meadat_t0.dat'; % data with time-zero correction

if ~exist(datafile,'file')
    error('You must go to the folder "pre-processed-air-2".'); 
end

% load a data set: 
data = dlmread(datafile); 

WaveSpeed = 0.3; % speed of light in free space m/ns. 
PropagationDist = 0.78; % distance in meters.
Nx = 51; Ny = 51; 

MaxTime = 10/0.3*WaveSpeed; % 10 nanoseconds, wave speed = 0.3 meters/ns. This is to make the wave speed in free space to be 1. 
[Nt,Np_mea] = size(data); 

t = linspace(0,MaxTime,Nt); % time variable
Nfreq = round((Nt+1)/2 + 0.25); % number of meaningful frequencies from Fourier transform 

X = linspace(-0.5,0.5,Nx)'; % grid in x
Y = X; % grid in y


% propagation in time domain: 
dataPropTime = wave_propagation_3d(reshape(data,Nt,Nx,Ny),PropagationDist,t,X,Y,WaveSpeed);
dataPropTime = reshape(dataPropTime,Nt,Np_mea);


% convert propagated time-domain data to frequency domain: 
MinFreq = 0; % minimum frequency
[dataPropFreq,freq] = FourierTransformFFT(dataPropTime,t,MinFreq,0,1);
dataPropFreq = dataPropFreq(1:Nfreq,:); % take only positive frequencies
freq = freq(1:Nfreq);



% ====================
% 2. propagate date in frequency domain. First, convert the original time
% domain data to frequency domain data, then propagate it. 


dataFreq = FourierTransformFFT(data,t,MinFreq,0,1);
dataFreq = dataFreq(1:Nfreq,:);

dataFreqProp = wave_propagation_freq_3d(reshape(dataFreq,Nfreq,Nx,Ny),PropagationDist,freq,X,Y,WaveSpeed); % data propagated in frequency domain
dataFreqProp = reshape(dataFreqProp,Nfreq,Np_mea)/2; % 1/2???



disp('Press any button to display the results:'); 
pause;
% display the results:
Nf = 120; % display only 20 first frequencies.

% MinReal = min(min(min(min(real(dataPropFreq)))), min(min(min(real(dataFreqProp)))));
% MinIm = min(min(min(min(imag(dataPropFreq)))), min(min(min(imag(dataFreqProp)))));
% MaxReal = max(max(max(max(real(dataPropFreq)))), max(max(max(real(dataFreqProp)))));
% MaxIm = max(max(max(max(imag(dataPropFreq)))), max(max(max(imag(dataFreqProp)))));

% figure; 
% for n = 1:Nf
%     subplot(1,2,1); 
%     surfc(real(reshape(dataPropFreq(n,:,:),Nx,Ny)')); 
%     title(sprintf('%s%5.2f%s','Frequency = ',freq(n),' GHz'));
%     axis([0 51 0 51 min(min(min(min(real(dataPropFreq(n,:,:))))),-200) max(max(max(max(real(dataPropFreq(n,:,:))))),200)]);
%     
%     view(70,15);
% 
%     colorbar; 
%     
%     subplot(1,2,2); 
%     surfc(real(reshape(dataFreqProp(n,:,:),Nx,Ny)')); 
%     title(sprintf('%s%5.2f%s','Frequency = ',freq(n),' GHz'));
%     axis([0 51 0 51 min(min(min(min(real(dataFreqProp(n,:,:))))),-200) max(max(max(max(real(dataFreqProp(n,:,:))))),200)]);
%     view(70,15);
%     colorbar; 
%     pause(0.2); 
% end
% 
% 
% figure; 
% for n = 1:Nf
%     subplot(1,2,1); 
%     surfc(imag(reshape(dataPropFreq(n,:,:),Nx,Ny)')); 
%     title(sprintf('%s%5.2f%s','Frequency = ',freq(n),' GHz'));
%     axis([0 51 0 51 min(min(min(min(imag(dataPropFreq(n,:,:))))),-200) max(max(max(max(imag(dataPropFreq(n,:,:))))),200)]);
%     
%     view(70,15);
% 
%     colorbar; 
%     
%     subplot(1,2,2); 
%     surfc(imag(reshape(dataFreqProp(n,:,:),Nx,Ny)')); 
%     title(sprintf('%s%5.2f%s','Frequency = ',freq(n),' GHz'));
%     axis([0 51 0 51 min(min(min(min(imag(dataFreqProp(n,:,:))))),-200) max(max(max(max(imag(dataFreqProp(n,:,:))))),200)]);
%     view(70,15);
%     colorbar; 
%     pause(0.2); 
% end
% 

n = 76; % index of frequency of 7.5GHz
figname = datafile(1:end-13);

figure; set(gca,'fontsize',15); 

imagesc(real(reshape(dataPropFreq(n,:,:),Nx,Ny)')); 
title(sprintf('%s%5.2f%s','Real part. Frequency = ',freq(n),' GHz'));
colorbar; 
xlabel('x'); ylabel('y'); 
print('-depsc2',[figname,'dataprop_time_re.eps']);

imagesc(imag(reshape(dataPropFreq(n,:,:),Nx,Ny)')); 
title(sprintf('%s%5.2f%s','Imaginary part. Frequency = ',freq(n),' GHz'));
colorbar; 
xlabel('x'); ylabel('y'); 
print('-depsc2',[figname,'dataprop_time_im.eps']);

imagesc(real(reshape(dataFreqProp(n,:,:),Nx,Ny)')); 
title(sprintf('%s%5.2f%s','Real part. Frequency = ',freq(n),' GHz'));
colorbar; 
xlabel('x'); ylabel('y'); 
print('-depsc2',[figname,'dataprop_freq_re.eps']);

imagesc(imag(reshape(dataFreqProp(n,:,:),Nx,Ny)')); 
title(sprintf('%s%5.2f%s','Imaginary part. Frequency = ',freq(n),' GHz'));
colorbar; 
xlabel('x'); ylabel('y'); 
print('-depsc2',[figname,'dataprop_freq_im.eps']);




% figure; 
% for n = 1:Nf
%     subplot(1,2,1); 
%     imagesc(real(reshape(dataPropFreq(n,:,:),Nx,Ny)')); 
%     title(sprintf('%s%5.2f%s','Real part. Data propagated in time domain. Frequency = ',freq(n),' GHz'));
%     colorbar; 
%     
%     subplot(1,2,2); 
%     imagesc(real(reshape(dataFreqProp(n,:,:),Nx,Ny)')); 
%     title(sprintf('%s%5.2f%s','Real part. Data propagated in frequency domain. Frequency = ',freq(n),' GHz'));
%     colorbar; 
%     pause(0.2); 
% end
% 
% figure; 
% for n = 1:Nf
%     subplot(1,2,1); 
%     imagesc(imag(reshape(dataPropFreq(n,:,:),Nx,Ny)')); 
%     title(sprintf('%s%5.2f%s','Imaginary part. Data propagated in time domain. Frequency = ',freq(n),' GHz'));
%     colorbar; 
%     
%     subplot(1,2,2); 
%     imagesc(imag(reshape(dataFreqProp(n,:,:),Nx,Ny)')); 
%     title(sprintf('%s%5.2f%s','Imaginary part. Data propagated in frequency domain. Frequency = ',freq(n),' GHz'));
%     colorbar; 
%     pause(0.2); 
% end
