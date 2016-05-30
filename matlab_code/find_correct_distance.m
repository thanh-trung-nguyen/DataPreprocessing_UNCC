function sol = find_correct_distance(dTR,dmea,a)

% find the correct distance from Tx to object when a signal not at the center is chosen
% using the Newton method
% the objective function is convex and has only 1 solution in
% [0,(dmea-dTR)/2].

% dOLD = (dmea - dTR)/2
% STOP = 0; Tolerance = 0.0001;
% 
% while ~STOP
%     fOLD = 2*dOLD*dTR + 2*dmea*sqrt(a^2 + dOLD^2) - dmea^2 + dTR^2;
%     derfOLD = 2*dTR + 2*dmea*dOLD/sqrt(a^2 + dOLD^2);
%     dNEW = dOLD - fOLD/derfOLD;
%     if abs(dNEW - dOLD) < Tolerance
%         STOP = 1;
%     end
%     dOLD = dNEW;
% end
% sol = dOLD;

%sol = (dmea^2 - dTR^2 - a^2)/2/(dmea+dTR);

if dTR < 0
    error('invalid input argument: dTR must be nonnegative');
end
d0 = 0; d1 = (dmea - dTR)/2; 
Tolerance = 0.0001;
while d1 - d0 > Tolerance
    d = (d0+d1)/2; 
    f = (d*sqrt((d+dTR)^2 + a^2) + sqrt((d+dTR)^4 + a^2*dTR^2))/(d+dTR) - dmea;
    if f < 0; 
        d0 = d; 
    else
        d1 = d;
    end  

end
sol = (d1 + d0)/2;