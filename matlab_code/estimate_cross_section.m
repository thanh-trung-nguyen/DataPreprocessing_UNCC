function MaskFile = estimate_cross_section(u,threshold)
% estimate the cross section from a data or a coefficient
% NOTE: the estimated cross section corresponds to the maximum values of u.
% Therefore if u is convex (has a local minimum), we need to change its sign in the input
% parameter. 
% threshold: a threshold value in percentage. 
% Nguyen Trung Thanh, 2014. 




% truncate the coefficient:
[Nx,Ny,Nz] = size(u);
max_value = max(max(max(u))); 
min_value = min(min(min(u)));
value_range = max_value - min_value;

MaskFile = zeros(Nx,Ny);

% % use the maximum layer: 
for i = 1:Nz
    M = max(max(u(:,:,i)));
    if M==max_value
        iz = i;
    end
end
u = u(:,:,iz);

MaskFile(u >= threshold*value_range + min_value) = 1;

if max(max(MaskFile)) == 0
    error('The cross-section is not shown up, check the input function');
end







