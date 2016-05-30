function data2 = offset_correction(data)
% function "offset_correction":
% does the offset correction by subtracting the mean value of the time dependent data. 
% INPUT: a 2d matrix of size MxN, M: number of samples in time, N: number of detectors.
% Output: a data of the same size as the input one.
% 
% @Nguyen Trung Thanh, UNCC 2013.
% NOTE: The function does not check for the correctness of the input data. 


Np = size(data,2);
for n = 1:Np
    dat = data(:,n); 
    data(:,n) = dat - mean(dat);
end
data2 = data;
