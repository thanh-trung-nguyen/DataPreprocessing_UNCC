function us = laplacetransform(ut,t,s)
% Calculate the laplace transform 
% size of ut: each column is a time signal. 
% t: can be a vector of time variable or just the time step. If it is the
% time step, ut should start from ut(dt,x), i.e., data at t =0 should not
% be in the input ut.
% 



Nt = size(ut,1);

if size(s,2) > 1
    s = s'; % s must be a column vector
end
if size(t,1) > 1
    t = t'; % t must be a row vector
end
if length(t) == 1
    dt = t; 
    t = linspace(dt,Nt*dt,Nt);
end
    
us = exp(-s*t)*ut*(t(2)-t(1));

