function dist = distance_estimation(u,mindist,dt,NoiseLevel,TxPos,RxPos)
% estimate the distance from the target to the receiver:
% u: a 1d time dependent curve
% mindist: minimum distance (an estimate) from the receiver to the target.
% dt: time step
% NoiseLevel: the data noise level
% TxPos, RxPos: positions of Tx and Rx in the down-range direction
% Nguyen Trung Thanh, UNCC 2013. 



MinIdx = round(mindist/dt) + 1;

u(1:MinIdx-1) = 0; 
u2 = abs(u);
u2(u2 < NoiseLevel) = 0;


idx = 1; % find the first index of nonzero of u2:
while (idx < length(u)-1) && (u2(idx) < eps)
    idx = idx + 1;
end
if idx == length(u)
    error('The input signal is too weak, choose another curve');
end
idx2 = idx-1; % find the first nonzero index of the first peak: 
while (u(idx)*u(idx2) > 0)
    idx2 = idx2 - 1;
end

dist = idx2*dt;
% dTR = RxPos - TxPos;
if nargin > 4
    dist = (dist + RxPos-TxPos)/2;
end

