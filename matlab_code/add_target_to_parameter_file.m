function  add_target_to_parameter_file(input_parfile,ObjType,ObjIdx,ObjCoord)
% this function add the location of the target inside the parameter files
% for forward solver of Larisa. 
% input_parfile: the parameter file that we want to add the location of target
% in.


output_parfile = [input_parfile(1:end-4),'_obj',num2str(ObjIdx),'.dat'];

fid1 = fopen(input_parfile,'r');
fid2 = fopen(output_parfile,'w');

% copy the first line:
tline = fgetl(fid1); 
fprintf(fid2,'%s\n',tline);

% copy the second and the third lines, modify the grid file name:

for i = 2:3
    text = fscanf(fid1,'%s',1);
    gridfilename = fscanf(fid1,'%s',1);
    newgridfilename = [gridfilename(1:end-4),num2str(ObjIdx),'.inp'];
    fprintf(fid2,'%s%s%s\n',text,'                        ',newgridfilename);
end
fgetl(fid1); 

% copy the input parameter file to the output parameter file:
for i = 4:23
    tline = fgetl(fid1); 
    %disp(tline); 
    fprintf(fid2,'%s\n',tline);
end; 

% update the target type: 
    text = fscanf(fid1,'%s',1);
    fprintf(fid2,'%s%s%s\n',text,'                        ',ObjType);
    fgetl(fid1);
    
for i = 25:26
    tline = fgetl(fid1); 
    %disp(tline); 
    fprintf(fid2,'%s\n',tline);
end; 


% update the coordinate: 
text = fscanf(fid1,'%s',1); %disp(text); 
fprintf(fid2,'%s%s%5.3f%s%5.3f%s%5.3f%s%5.3f%s%5.3f%s%5.3f\n',text, '        ', ObjCoord(1),'  ', ObjCoord(2),'  ', ObjCoord(3),'  ', ObjCoord(4),'  ', ObjCoord(5),'  ', ObjCoord(6));
    
fclose(fid1);
fclose(fid2);
