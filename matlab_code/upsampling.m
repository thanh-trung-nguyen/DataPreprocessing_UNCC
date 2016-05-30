function [intdat,upstime] = upsampling(data,propagation_vector,timefactor,method)
% interpolate the data w.r.t. time domain
% method used: spline or zero padding
% Syntax: 
%     intdat = upsampling(data,propagation_vector,timefactor)
%     intdat = upsampling(data,propagation_vector,timefactor,method)
%     [intdat,upstime] = upsampling(data,propagation_vector,timefactor,method)
% 
% INPUT: 
%     data: A-scan or B-scan or multiple B-scans
%     propagation_vector: acquisition time
%     timefactor: factor for upsampling (e.g. 4 times)
%     method: the method used, must be either 'spline' or 'zeropadding'.
%     
% OUTPUT:
%     intdat: interpolated data
%     upstime: upsampled propagation vector
%     
% *************************************************************************
% *************************************************************************

% check the input data: 
if nargin < 3
    error('Not enough input arguments');
end
if size(data,1) ~= length(propagation_vector)
    error('The data and the propagation vector do not match!');
end
if nargin < 4
    method = 'zeropadding'; % the default method
end


[Nt,Nx,Nr] = size(data);
upstime = propagation_vector(1):(propagation_vector(end)-propagation_vector(1))/((Nt-1)*timefactor):propagation_vector(end);
NtNew = length(upstime);
intdat = zeros(NtNew,Nx,Nr);
if strcmp(method,'spline')
    for k = 1:Nr
        intdat(:,:,k) = (spline(propagation_vector,data(:,:,k)',upstime))';
    end
elseif strcmp(method,'zeropadding')       
    for k = 1:Nr
        for j = 1:Nx
            ascan = data(:,j,k);
            fa = fft(ascan);
            ExtAscan = zeros(timefactor*(Nt-1)+1,1);
            halflength = floor(Nt/2);
            ExtAscan(1:halflength) = fa(1:halflength);
            ExtAscan(end-halflength:end) = fa(end-halflength:end);
            intdat(:,j,k) = real(ifft(ExtAscan)*timefactor);
        end
    end    
elseif strcmp(method,'linear')    
    for k = 1:Nr    
        i = 1:Nt;        
        intdat((i-1)*timefactor+1,:,k) = data(i,:,k);
        i = i(1:end-1);
        for j = 2:timefactor
            weight = (j-1)/timefactor; 
            intdat((i-1)*timefactor+j,:,k) = data(i,:,k)*(1 - weight) + data(i+1,:,k)*weight;
        end
    end    
else
    error('The method must be either linear or zeropadding');
end
