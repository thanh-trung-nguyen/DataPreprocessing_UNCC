function [threshold,threshold_percent] = estimate_truncation_value(C,Size)
% estimate the truncation threshold in the coefficient value based on the
% given xy cross-section size. The size is counted from the maximum value
% of C. 

% C: a 3D or 2D data set
% size: a given size, with the same unit as the data. 


max_value = max(max(max(C))); 
min_value = min(min(min(C)));
value_range = max_value - min_value; 

dt = 0.01; Threshold = dt:dt:1-dt;
N = length(Threshold);

s = zeros(N,1);

for i = 1:N
    C1 = C; 

    C1(C < min_value + Threshold(i)*value_range) = 0;
    C1(C >= min_value + Threshold(i)*value_range) = 1;

    C1 = sum(C1,3);
    C1(C1 > 0) = 1;
    s(i) = sum(sum(C1));
end
s = abs(s - Size);

threshold_percent = mean(Threshold(s==min(s)));
threshold = min_value + threshold_percent*value_range;
