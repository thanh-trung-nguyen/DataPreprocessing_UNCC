function [datF,freq] = data_time2freq(datfile, t, Nx, Ny,mainFreq)
% Convert from time domain radar data to frequency domain data using
% Fourier transform. 
% Note: the frequencies of the transformed data are determined from the
% Fast Fourier transform, they are not arbitrary. 
% Input: datfile: the data file in time domain, each column is the
% time-dependent data at a certain position of radar.
% t: a vector of time steps of the data. 
% Optional: Nx, Ny: number of measurement grid points in x and y
% directions, this is only for displaying the data
% mainFreq: the central frequency of the data, in GHz. 
% Output: 
% datF: data in frequency domain. Each column is the data at a certain radar position. 
% freq: a vector of frequencies. 
% 
% @Thanh Nguyen, Iowa State University, 2016. 



dat = dlmread(datfile);
[Nt,Npos] = size(dat); % Nt: number of samples in time, Nx: number of radar positions.

datF = zeros(Nt,Npos);  % data in frequency domain
for n = 1:Npos
    [datF(:,n),freq]  = fourier_fft(dat(:,n),t); 
end

freq = freq(ceil(Nt/2):end); % just take positive frequencies, the data is symmetric. 
datF = datF(ceil(Nt/2):end,:); 

if nargin > 2
    
    % ploting the results: optional
    figname = datfile(1:end-4); %  save the figures in the same folder as data, but with different names
    f = abs(freq - mainFreq);
    Idx = find(f==min(f)); % find the index of the main frequency


    for n = -10:10 % print 3 figures at 3 frequencies
        CurrentIdx = Idx + n*5; 
        currentFreq = freq(CurrentIdx);
        figure(1); 
        set(gca,'fontsize',15); 
        imagesc(real(reshape(datF(CurrentIdx,:),Nx,Ny))); colorbar; 
        title(sprintf('%s%5.1f%s','Real part of the data. Frequency = ', currentFreq,' GHz')); 
        pause(1);

    end

    for n = -1:1 % print 3 figures at 3 frequencies
        CurrentIdx = Idx + n*20;
        currentFreq = freq(CurrentIdx);
        figure(1); 
        set(gca,'fontsize',15); 
        imagesc(real(reshape(datF(CurrentIdx,:),Nx,Ny))); colorbar; 
        title(sprintf('%s%5.1f%s','Real part of the data. Frequency = ', currentFreq,' GHz')); 
        xlabel('x'); ylabel('y'); 
        print('-depsc2',sprintf('%s%s%3.1f%s',figname,'_freq',currentFreq,'_real.eps'));

        imagesc(imag(reshape(datF(CurrentIdx,:),Nx,Ny))); colorbar; 
        title(sprintf('%s%5.1f%s','Imaginary part of the data. Frequency = ', currentFreq,' GHz')); 
        xlabel('x'); ylabel('y'); 
        print('-depsc2',sprintf('%s%s%3.1f%s',figname,'_freq',currentFreq,'_imag.eps'));
    end
end