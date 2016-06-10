% frequency data propagation for the sandbox


% % -------------load data of object and sand alone: change the path to your
% % actual data files:
DisplayedFreqList = [2.6e9, 4.7e9, 7e9]; % list 
% % datafile = '/home/thanhnguyen/Programming_Git/UNCC_data/FreqData_2016/non-blind/Sandbox/Data_files/sandbox-aluminum-cylinder.csv';
% % datafile = 'C:/Users/thanh/Thanh-data/work/Programming_Git/UNCC_data/FreqData_2016/non-blind/Sandbox/Data_files/sandbox-aluminum-cylinder.csv';
% datafile = 'C:/Users/thanh/Thanh-data/work/Programming_Git/UNCC_data/FreqData_2016/non-blind/Sandbox/Data_files/sandbox-water-bottle.csv';
% datafile = '/home/thanhnguyen/Programming_Git/UNCC_data/FreqData_2016/non-blind/Sandbox/Data_files/sandbox-sand-only.csv';
% [freq,data3d,data2d] = load_data_freq(datafile,300,50,51);

figname = 'sand-only'; % name of figures to save on computer

% datafile = 'sandbox-aluminum-cylinder.csv';

DisplayedFreqList = [2.6e9, 4.7e9, 7e9]; % list 


DisplayedFreqList = [2.6e9, 4.7e9, 7e9]; % li
DisplayedFreqList = [2.6e9, 4.7e9, 7e9]; % list st of frequency to visualize the result: 

d = 0.75; % distance of propagation. propagate to the front surface of the sand box.

% -----------the above parameters may be changed. 
% -----------the above parameters may be changed. 

NfreqDisplayed = length(DisplayedFreqList);
[Nfreq,Nx,Ny] = size(data3d);
% propagate data:
LightSpeed = 299792458; % m/s.
X = linspace(-0.5,0.5,Nx); 
Y = linspace(-0.5, 0.5, Ny); 


% data propagation in frequency domain: 
dataProp = wave_propagation_freq_3d(data3d,d,freq',X,Y,LightSpeed); 

dataProp = reshape(dataProp,Nfreq,Nx*Ny);

% figure;

% -----------the above parameters may be changed. 
DisplayedFreqList = [2.6e9, 4.7e9, 7e9]; % list % for n = 1:Nfreq; 
%     image
DisplayedFreqList = [2.6e9, 4.7e9, 7e9]; % list sc(abs(reshape(dataProp(n,:),50,51))'); 
%     title(freq(n)); colorbar;  
%     pause(0.2); 
% end

imagesc(real(data2d));
set(gca,'ytick',[1:20:Nfreq]);
set(gca,'yticklabel',freq(1:20:end)/1e9);
xlabel('Measurement point');
ylabel('Frequency (GHz)');
print('-depsc2',sprintf('%s%s',figname,'origdat_re.eps'));

imagesc(real(dataProp));
set(gca,'ytick',[1:20:Nfreq]);
set(gca,'yticklabel',freq(1:20:end)/1e9);
xlabel('Measurement point');
ylabel('Frequency (GHz)');
print('-depsc2',sprintf('%s%s%4.2f%s',figname,'_dist',d,'_re.eps'));


% save some figures: 
for k = 1:NfreqDisplayed

    x = abs(freq - DisplayedFreqList(k));
    n = find(x ==min(x),1,'first');

    % the original data: 
    imagesc(real(reshape(data3d(n,:,:),Nx,Ny))'); 
    title(sprintf('%s%5.2f%s','Real part. Frequency = ',freq(n)/10^9,' GHz'));
    colorbar; 
    xlabel('x'); ylabel('y'); 
    print('-depsc2',sprintf('%s%s%4.2f%s',figname,'_origdat_freq',freq(n)/10^9,'_re.eps'));

    % propagated data:
    imagesc(real(reshape(dataProp(n,:,:),Nx,Ny))'); 
    title(sprintf('%s%5.2f%s%4.2f%s','Real part. Frequency = ',freq(n)/10^9,' GHz. Propagation distance = ',d,' (m)'));
    colorbar; 
    xlabel('x'); ylabel('y'); 
    print('-depsc2',sprintf('%s%s%4.2f%s%4.2f%s',figname,'_dist',d','_freq',freq(n)/10^9,'_re.eps'));

    % the original data: 
    imagesc(abs(reshape(data3d(n,:,:),Nx,Ny))'); 
    title(sprintf('%s%5.2f%s','Modulus. Frequency = ',freq(n)/10^9,' GHz'));
    colorbar; 
    xlabel('x'); ylabel('y'); 
    print('-depsc2',sprintf('%s%s%4.2f%s',figname,'_origdat_freq',freq(n)/10^9,'_abs.eps'));

    % propagated data:
    imagesc(abs(reshape(dataProp(n,:,:),Nx,Ny))'); 
    title(sprintf('%s%5.2f%s%4.2f%s','Modulus. Frequency = ',freq(n)/10^9,' GHz. Propagation distance = ',d,' (m)'));
    colorbar; 
    xlabel('x'); ylabel('y'); 
    print('-depsc2',sprintf('%s%s%4.2f%s%4.2f%s',figname,'_dist',d','_freq',freq(n)/10^9,'_abs.eps'));


end
