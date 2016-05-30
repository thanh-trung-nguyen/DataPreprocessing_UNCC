function add_two_targets_to_parameter_file(parfile_common, ObjIdx, inpfile1, inpfile2)

outputfile = [parfile_common(1:end-4),'_obj',num2str(ObjIdx),'.dat'];

fid1 = fopen(inpfile1,'r');
fid2 = fopen(inpfile2,'r');
fid0 = fopen(parfile_common,'r');

fidout = fopen(outputfile,'w');

% copy the first line:
tline = fgetl(fid0); 
fprintf(fidout,'%s\n',tline);

% copy the second and the third lines, modify the grid file name:

for i = 2:3
    text = fscanf(fid0,'%s',1);
    gridfilename = fscanf(fid0,'%s',1);
    newgridfilename = [gridfilename(1:end-4),num2str(ObjIdx),'.inp'];
    fprintf(fidout,'%s%s%s\n',text,'                        ',newgridfilename);
end
fgetl(fid0); 

% copy the input parameter file to the output parameter file:
for i = 4:22
    tline = fgetl(fid0); 
    %disp(tline); 
    fprintf(fidout,'%s\n',tline);
end; 

% update the number of targets: 
    text = fscanf(fid0,'%s',1);
    fprintf(fidout,'%s%s%d\n',text,'                        ',2);
    fgetl(fid0);

% go to the target's properties of the two halves:
for i = 1:23
    fgetl(fid1);
    fgetl(fid2);
end

% copy the target's properties from target 1:
for i = 24:27
    tline = fgetl(fid1); 
    %disp(tline); 
    fprintf(fidout,'%s\n',tline);
end; 

% copy the target's properties from target 2:
for i = 24:27
    tline = fgetl(fid2); 
    %disp(tline); 
    fprintf(fidout,'%s\n',tline);
end; 
    
fclose(fid0);
fclose(fidout);    
fclose(fid1);
fclose(fid2);




