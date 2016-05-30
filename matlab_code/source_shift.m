function data2 = source_shift(data,distTxObj_new,distRxObj,dt,distTxObj_old)
%
% Move the source in such a way that the new data corresponds to a given distance for inversion. 
% % @Nguyen Trung Thanh, 2013.
% -----------

Threshold = 0.8;
[Nt,Np] = size(data); 

if nargin<5
    % estimate the current time delay: 
    [idt,idx] = find(data==min(min(data))); %#ok<ASGLU>
    idx = idx(1); 
    u = data(:,idx);

    idt = find_first_sample(u,Threshold); % find the first sample of the target's signal at the strongest detector

    % calculate the new time delay: 
    newIdt = round((distTxObj_new + distRxObj)/dt) + 1;

    NrSampleShift = idt - newIdt
else
     NrSampleShift = round((distTxObj_old - distTxObj_new)/dt)
end

data2 = zeros(Nt,Np);

if NrSampleShift >= 0
     data2(1:Nt-NrSampleShift,:) = data(NrSampleShift+1:Nt,:);
else
     data2(NrSampleShift+1:Nt,:) = data(1:Nt-NrSampleShift,:);
end