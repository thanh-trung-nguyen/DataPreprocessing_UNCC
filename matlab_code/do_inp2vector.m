function do_inp2vector(ObjIdx,inputfileprefix,outputfileprefix)
% run the function inp2vector.m for different data sets at once
% example: do_inp2vector(1:10,'object','object');


for i = 1:length(ObjIdx)
    inputfile = [inputfileprefix,num2str(ObjIdx(i)),'.inp'];
    outputfile = [outputfileprefix,num2str(ObjIdx(i)),'.m'];
    if exist(inputfile,'file')
        inp2vector(inputfile,outputfile);
    else
        fprintf('%s%s\n','File: ', inputfile,' does not exist');
    end
end
    


