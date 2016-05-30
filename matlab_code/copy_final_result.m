function copy_final_result(ObjIdx,folder,final_result_folder,objfilename)
% copy the final results of all targets into the same folder for
% postprocessing:
% For each data set (object), there will be two files with extensions .dat
% and .inp. The first file contains the coefficient values at the FEM
% nodes. The second file is the coefficient in the structure of FEM grid.
% we can visualize the reconstruction results from this folder. 

if nargin < 4
    objfilename = 'object';
end
if nargin < 3
    final_result_folder = 'final_results';
end


Nobj = length(ObjIdx);
if nargin > 1
    eval(['cd ' folder]);
end
if ~exist(final_result_folder,'dir')
    eval(['mkdir ' final_result_folder]);
end

for k = 1:Nobj
        i = ObjIdx(k);
        
        fname1 = ['object_',num2str(i),'/Test_c.m'];
        if exist(fname1,'file')
            
%             fname2 = ['final_results/object_',num2str(i),'.dat'];

            err = dlmread(['object_',num2str(i),'/Test_L2Error.m']);              
            [nx,ny] = find(err == min(min(err)));
            nx = nx(1);
            ny = ny(1);

%             err = reshape(err',1,size(err,1)*size(err,2));
%             idx = find(err == min(err),1,'first');
%             Coeff = dlmread(fname1);
%             Coeff = Coeff(idx,:);
%             dlmwrite(fname2,Coeff,'delimiter',' ');


            fname3 = ['object_',num2str(i),'/Coef_',num2str(nx),'_',num2str(ny-1),'.inp'];
            fname4 = [final_result_folder,'/',objfilename,num2str(i),'.inp'];
            [~,message] = copyfile(fname3,fname4);
            disp(message);
        end        
end