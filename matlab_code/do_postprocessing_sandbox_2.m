ObjIdx = 3:8;
datafolder ='/home/tnguy152/UNCC_2012/programs/WORK_GCA_sandbox_2';
datafile = 'prop_dat_inv_w_s9.00_7.00_41.m';
final_result_folder = 'final_results_Larisa';
% final_result_folder = 'final_results';
objfilename = 'object_'; % common part of object result files to be created


parfile_fp = 'parameter_forprob.dat';
parfile_vis = 'parameter_vis.dat'; % parameter file for the visualization
% parfile_preprocessing = '/home/thanhnguyen/UNCC_2012/programs/data_preprocessing/examples/data_sandbox/parameters_for_preprocessing.dat';
parfile_preprocessing = '/home/tnguy152/UNCC_2012/programs/WORK_GCA_sandbox_2/parameters_for_preprocessing.dat';

bckgr_coef = 4; % back ground dielectric constant

Threshold_CrossSection = 0.4; % threshold in estimating the cross section from data. 
% This value should be provided based on the estimation of the cross section of the calibrating target

thickness_threshold_Larisa = 0.9; % threshold in estimating the thickness of targets used in Larisa's tests.


eval(['cd ' datafolder]);

% %% step 1: copy the final result to a folder named "final_results":
% copy_final_result(ObjIdx,datafolder,final_result_folder,objfilename);
% 
% %% step 2: extract the coefficient values from the INP file:
% INP_file_prefix = [final_result_folder,'/',objfilename];
% do_inp2vector(ObjIdx,INP_file_prefix, INP_file_prefix);        

%% step 3: estimate the cross-section, truncate the reconstructed coefficient and create the parameter
%% files for the visualization. 

truncate_reconstructed_coefficient(bckgr_coef,datafolder,datafile,final_result_folder,objfilename,parfile_preprocessing,...
    parfile_fp,parfile_vis,object_properties_sandbox_2,...
    ObjIdx,Threshold_CrossSection,thickness_threshold_Larisa,targetdepths_sandbox_2);


