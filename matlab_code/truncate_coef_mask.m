function newdata = truncate_coef_mask(coef,mask)
% truncate a 3d coefficient in the first and second directions based on the
% mask


newdata = coef*0 + 1; 

for iz = 1:size(coef,3)
    u = coef(:,:,iz);
    u(mask==0) = 1;
    newdata(:,:,iz) = u;
end