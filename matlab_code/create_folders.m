function create_folders(mainfolder,subfolder,StartIdx,EndIdx)
% Function create_folders.m
%
% create a folder named mainfolder, and subfolders with prefix subfolder
% and numbered from StartIdx to EndIdx
% for example: to create the main folder "data_set_1" in the current directory. The subfolders are
% "object_1", "object_2", "object_3", we call
% create_folders('data_set_1','object_',1,3);
% Nguyen Trung Thanh, UNCC 2014
%
% =========================================================================


eval(['mkdir ' mainfolder]);
eval(['cd ' mainfolder]);
if StartIdx <= EndIdx
    for i = StartIdx:EndIdx
        eval(['mkdir ' [subfolder,num2str(i)]]);
    end
end
cd ..

    



