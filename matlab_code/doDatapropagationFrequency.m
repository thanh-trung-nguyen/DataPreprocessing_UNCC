function doDatapropagationFrequency(datafolder,datafilename,incwave,PropDistance,DisplayedFreqList,X,Y)
% run data propagation of the frequency domain data. The main computation
% is done in the routine "wave_propagation_freq_3d.m".
% 
% INPUT: 
%     1. datafolder: the folder containing the data file to be propagated.
%     It can be relative or absolute path.
%     2. datafilename: the name of the data file, without folder name
%     3. incwave: a 3D matrix of incident wave (data without any target). 
%     The first index is in frequency, the second index is in x, the third index is in y. 
%     This incident wave matrix can be loaded using load_data_freq.m function. 
%     4. PropDistance: propagation distance in meters, could be a vector. 
%     5. DisplayedFreqList: a list of frequencies at which the data is displayed.
%     Note that this is not the list of all frequencies in the measured
%     data.
%     6, 7: X, Y (optional): the grid in the x and y directions. If not
%     provided, the default vectors of X = linspace(-0.5,0.5, Nx), Y =
%     linspace(-0.5,0.5,Ny) are used, where Nx, Ny are the number of grid
%     points in x and y directions given by the incident wave matrix.
%
% OUTPUT:
%     This function will create figures of real, imaginary components, and the modulus of the propagated data. 
%     The output figures are in EPS format for embedding into latex
%     reports. The figures are saved in the current folder. Their names
%     start with the name of the data file.
%
% Modified on May 30, 2016 by Thanh Nguyen. 
% ============================================
% The following parameters can be used: 
% FreqList = [2.09e9, 2.45e9, 2.91e9, 4.55e9, 7e9]; % choose indices of frequencies to be displayed for the above distance of propagation
% /home/thanhnguyen/Programming_Git/UNCC_data/FreqData_2016/non-blind/Objects-in-air/Data-files/



LightSpeed = 299792458; % m/s.


datafile = [datafolder, datafilename];
fprintf(['Data file: ',datafile,'\n']);
[Nfreq,Nx,Ny] = size(incwave); 
NfreqDisplayed = length(DisplayedFreqList);

% load the data into Matlab:  
if ~exist(datafile,'file');
    error('Data file is not found. Please check the path');
end
[freq,data3d] = load_data_freq(datafile,Nfreq,Nx,Ny);
% note: freq is in Hz. 

% extract the scattered wave by subtracting the incident wave from the
% total wave:
data3d = data3d - incwave;

% ------the default grids in X and Y: these are for UNCC frequency data.
if nargin < 6    
    Xmin = -0.5; Xmax = 0.5; 
    Ymin = -0.5; Ymax = 0.5; 
    X = linspace(Xmin,Xmax,Nx); 
    Y = linspace(Ymin,Ymax,Ny);
end


% ---- propagate data and dispay results:
for id = 1:length(PropDistance)
    d = PropDistance(id); % current distance of propagation
    
    dataProp = wave_propagation_freq_3d(data3d,d,freq',X,Y,LightSpeed);
    dataProp = reshape(dataProp,Nfreq,Nx*Ny);


    % --- display the result and save some figures:
    
    
    idx = strfind(datafilename,'.'); idx = idx(end); % find the dot to remove the extension. 
    figname = datafilename(1:idx-1);

    
    figure; set(gca,'fontsize',15); 
    
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
        surf(abs(reshape(dataProp(n,:,:),Nx,Ny))'); 
        title(sprintf('%s%5.2f%s%4.2f%s','Modulus. Frequency = ',freq(n)/10^9,' GHz. Propagation distance = ',d,' (m)'));
        colorbar; 
        xlabel('x'); ylabel('y'); 
        print('-depsc2',sprintf('%s%s%4.2f%s%4.2f%s',figname,'_dist',d','_freq',freq(n)/10^9,'_abs.eps'));


    end
end

