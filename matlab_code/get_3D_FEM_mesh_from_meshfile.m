function [Xmesh,Ymesh,Zmesh] = get_3D_FEM_mesh_from_meshfile(inputfile)
% this function get the coordinates of a 3D FEM mesh from a given mesh file
% in the format INP. 
% 
% Nguyen Trung Thanh, 2014.


fid = fopen(inputfile);
if fid == -1
    error('The input file does not exist');
end

nr_nodes = round(str2double(fscanf(fid,'%s',1)));
fgetl(fid);


Xmesh = zeros(1,nr_nodes); Ymesh = Xmesh; Zmesh = Xmesh;
for i =1:nr_nodes
    idx = round(str2double(fscanf(fid,'%s',1)));
    Xmesh(idx) = str2double(fscanf(fid,'%s',1));
    Ymesh(idx) = str2double(fscanf(fid,'%s',1));
    Zmesh(idx) = str2double(fscanf(fid,'%s',1));    
end


fclose(fid);

