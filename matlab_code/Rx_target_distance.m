function dist = Rx_target_distance(u,mindist,dt,NoiseLevel,TxPos,RxPos)
% estimate the distance from the target to the receiver:
% u: a 1d time dependent curve
% mindist: minimum distance (an estimate) from the receiver to the target.
% dt: time step
% NoiseLevel: the data noise level
% TxPos, RxPos: positions of Tx and Rx in the down-range direction
% Nguyen Trung Thanh, UNCC 2013. 



MinIdx = round(mindist/dt) + 1;

u(1:MinIdx-1) = 0; 
% find the first signal which is stronger than the noise level: 
idx = find(u < -NoiseLevel,1,'first');
if idx == length(u)
    error('The input signal is too weak, choose another curve');
end

% go back to the negative peak before this peak:
while idx > 1 && u(idx) <= 0
    idx = idx - 1;
end
while idx > 1 && u(idx) >0
    idx = idx - 1;
end
while idx > 1 && u(idx) <= 0
    idx = idx - 1;
end
idx = idx + 1;

dist = idx*dt;
% dTR = RxPos - TxPos;
if nargin > 4
    dist = (dist + RxPos-TxPos)/2;
end

