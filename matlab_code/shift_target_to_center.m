function [coef_new,masknew] = shift_target_to_center(coef_old,mask,bckgr_coef)

masknew = 0*mask;
coef_new = coef_old*0 + bckgr_coef; 

[Nx,Ny,Nz] =size(coef_old);

stat = regionprops(mask);

stat.Centroid(1)
stat.Centroid(2)
stat.BoundingBox(1)
stat.BoundingBox(2)


shift_column = round(Ny/2 - stat.Centroid(1));
shift_row = round(Nx/2 - stat.Centroid(2));

row_old = stat.floor(BoundingBox(2)): stat.floor(BoundingBox(2))+stat.BoundingBox(4);
col_old = stat.floor(BoundingBox(1)): stat.floor(BoundingBox(1))+stat.BoundingBox(3);

masknew(row_old+shift_row, col_old+shift_column) = mask(row_old,col_old);

for iz = 1:Nz
    coef_new(row_old+shift_row, col_old+shift_column,iz) = coef_old(row_old,col_old,iz);
end


