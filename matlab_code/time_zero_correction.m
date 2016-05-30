function data2 = time_zero_correction(data,distTxRx,dt,noiselevel,shift)
% function: time_zero 
% estimate the time zero and shift the data accordingly. 
% The direct signals from the transmitter to receiver is used as the reference
% Input: 
% data: a matrix of data. Each column is a time-domain curve. The data is measured in the xy plane
% distTxRx: distance between Tx and Rx
% dt: time step of the recorded signal 
% noiselevel: level of noise in the data, this is used to exclude data
% before the first peak of the direct signal. 
% shift: This parameter is used ocassionally only, in the case when the
% "first" direct signal, which goes backwards from the transmitter to the
% detector, is not detectable. In this case we choose the signal comming
% from the front of the horn to the detector as the first direct signal and
% then shift by a certain number of samples. This number is the same for
% all data sets. We have found that shift = 82 for the UNCC data (164 in the upsampled data). 
% 
% OUTPUT: data2: another matrix of the same size as the input data, with corrected time zero
% @Nguyen Trung Thanh, 2013
% 
% Remark: work for UNCC data with a fixed transmitter, not for bistatic data!
% Remark: this code assumes that the first peak is positive. To change to
% the first negative peak, comment line 42. 
% -----------------


[Nt,Np] = size(data); %Np: number of receiver positions
Ncenter = round((Np+1)/2); % take the central receiver
TrueTimeZero = round(distTxRx/dt) + 1;

if nargin > 4
    TrueTimeZero = TrueTimeZero + shift;
end

% take the curve at the center of the measurement plane (closest to the transmitter)
u = data(:,Ncenter); 
 
% find the first direct signal:
it = find(u < -noiselevel,1,'first'); % find first negative value with amplitude larger than the noise level
while u(it) < 0; it = it - 1;  end % move back till non-negative value
while u(it) >= 0; it = it - 1; end % move back till non-positive value. The incident wave has the first peak upwards!

TimeZero = it+1; % choose the first non negative value as the first sample of the direct signal
  
data2 = 0*data;
if TimeZero <= TrueTimeZero
    data2(TrueTimeZero-TimeZero+1:Nt,:)  = data(1:Nt-TrueTimeZero+TimeZero,:);    
else
    data2(1:Nt+TrueTimeZero-TimeZero,:)  = data(TimeZero-TrueTimeZero+1:Nt,:);        
end