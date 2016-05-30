function coef = inp2vector(inputfile,outputfile)
% this function extract the values of the coefficient in an INP file and
% write to a new file with the file name outputfile. 
% the input file is a file of finite element mesh. 


fid = fopen(inputfile);
if fid == -1
    error('The input file does not exist');
end

nr_nodes = round(str2double(fscanf(fid,'%s',1)));
nr_elements = round(str2double(fscanf(fid,'%s',1)));

fgetl(fid);

for i= 1:nr_nodes+nr_elements+2
    fgetl(fid);
end

coef = zeros(1,nr_nodes);
for i =1:nr_nodes
    fscanf(fid,'%s',1);
    coef(i) = str2double(fscanf(fid,'%s',1));
end

if nargin > 1
    dlmwrite(outputfile,coef,'delimiter',' ');
end



