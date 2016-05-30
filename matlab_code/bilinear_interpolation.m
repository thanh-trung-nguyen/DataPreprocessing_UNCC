function dataNew = bilinear_interpolation(data,Xold,Yold,Xnew,Ynew)
% function bilinear_interpolation:
% Compute the bilinear interpolation of a matrix.
% Nguyen Trung Thanh, UNCC 2013.
% *************************************************************************

EPS = eps;

N1 = length(Xold);
N2 = length(Yold);
N1new = length(Xnew);
N2new = length(Ynew);

% check the consistency of input data: 
if (N1 ~= size(data,1) || N2~= size(data,2))
    error('Input parameters are not consistent in length');
end
if (Xold(1) > Xnew(1) + EPS) || (Xold(N1) < Xnew(N1new) - EPS)
    error('Xnew is out of range'); 
end
if (Yold(1) > Ynew(1) + EPS) || (Yold(N2) < Ynew(N2new) - EPS)
    error('Ynew is out of range'); 
end


% interpolation: 
dataNew = zeros(N1new,N2new) + 1;
idx = 2;  % indices in the old grid

for i = 1:N1new
    while (idx < N1) && (Xnew(i) > Xold(idx)) %find the old index in x
        idx = idx + 1; 
    end
    idy = 2;
    for j = 1:N2new
        while (idy < N2) && (Ynew(j) > Yold(idy))
            idy = idy + 1;
        end        
        dataNew(i,j) = ( data(idx-1,idy-1)*(Xold(idx) - Xnew(i))*(Yold(idy) - Ynew(j)) ...
				        +data(idx-1,idy)*(Xold(idx) - Xnew(i))*(Ynew(j) - Yold(idy-1)) ...
				        +data(idx,idy-1)*(Xnew(i) - Xold(idx-1))*(Yold(idy) - Ynew(j)) ...
				        +data(idx,idy)*(Xnew(i) - Xold(idx-1))*(Ynew(j) - Yold(idy-1)) ...
				       )/((Xold(idx) - Xold(idx-1))*(Yold(idy) - Yold(idy-1)));			            
    end
end