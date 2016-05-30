function [freq,data3d,data2d] =  load_data_freq(datafile,NumFreq,Nx,Ny)
% LOAD_DATA_FREQ loads the frequency domain data. 
%Input: datafile: a CSV data file
%NumFreq: number of frequencies
%Nx: number of data points in the x (horizontal) direction
%Ny: number of data points in the y (vertical) direction
%
%Output: 
% freq: a row vector of frequencies
% data3d: a NumFreq x Nx x Ny matrix of the frequency domain data
% data2d: a 2D matrix of the data, data at each point is given by a column. 
% the row index is counted row by row in the original 3D data.
%
% @Thanh Nguyen, Iowa State University, 2016. 

% default number of frequencies, Nx, Ny. These numbers are given from the UNCC experimental setup. 
if nargin < 2
  NumFreq = 300;
  Nx = 50;
  Ny = 51; 
end

freq = zeros(1,NumFreq);
data3d = zeros(NumFreq,Nx,Ny);

fileID = fopen(datafile);
if ~fileID
  error('The file does not exist');
end

for iy = 1:Ny
  for ix = 1:Nx
    % load the three header lines:      
    line = fgets(fileID);
    line = fgets(fileID);
    line = fgets(fileID);
    
    for ifreq = 1:NumFreq
        string = fscanf(fileID,'%s',1);
        idx = find(string==';');
        
        freq(ifreq) = str2double(string(1:idx(1)-1));
        data3d(ifreq,ix,iy) = str2double(string(idx(1)+1:idx(2)-1)) ...
                        + 1i*str2double(string(idx(2)+1:idx(3)-1));
                    
      line = fgets(fileID);
      
    end  
  end
end
fclose(fileID);

% reshape the data so that it becomes a 2D matrix: 
if nargout > 2
%     data2d = zeros(NumFreq,Nx*Ny); 
%     for idy = 1:Ny
%         for idx = 1:Nx
%             col = idx + (idy-1)*Nx;
%             data2d(:,col) = data3d(:,idx,idy);
%         end
%     end
    data2d = reshape(data3d,NumFreq,Nx*Ny);

end

