function [x,y,z]=get_FEM_domain(parameter_file)
% get the FEM domain in x, y and z from a parameter file for the forward
% solver of waves24
% INPUT: parameter_file: a structured parameter file for the forward solver
% @Nguyen Trung Thanh, 2014


fid1 = fopen(parameter_file);

if ~fid1
    error('The input file cannot be opened');
end

% copy the input parameter file to the output parameter file:
for i = 1:5
    fgetl(fid1); 
end; 

fscanf(fid1,'%s',1); 

% load the z coordinates of the target in the input parameter files:
Nx = round(str2double(fscanf(fid1,'%s',1)));
Ny = round(str2double(fscanf(fid1,'%s',1)));
Nz = round(str2double(fscanf(fid1,'%s',1)));

dx = str2double(fscanf(fid1,'%s',1));
dy = str2double(fscanf(fid1,'%s',1));
dz = str2double(fscanf(fid1,'%s',1));

Xmin = str2double(fscanf(fid1,'%s',1));
Ymin = str2double(fscanf(fid1,'%s',1));
Zmin = str2double(fscanf(fid1,'%s',1));

x = (0:Nx-1)*dx+Xmin;
y = (0:Ny-1)*dy+Ymin;
z = (0:Nz-1)*dz+Zmin;

fclose(fid1);

