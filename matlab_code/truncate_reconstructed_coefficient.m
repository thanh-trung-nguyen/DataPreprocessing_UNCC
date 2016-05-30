function truncate_reconstructed_coefficient(bckgr_coef,datafolder,datafilename,final_result_folder,objfilename,parfile_preprocessing,parfile_forprob,parfile_vis,object_properties_func,ObjIdx,size_threshold,thickness_threshold,target_depth_file)
% estimate the cross-section and the location of the targets in the
% reconstructed coefficient. 
% also create a parameter file with the updated target's properties for the
% visualization. 
% INPUT: 
% bckgr_coef: the dielectric constant of the background medium
% datafolder: the main folder containing the object subfolders and the final_results folder
% datafilename: the data file in the Laplace transform
% parfile_preprocessing: file containing the parameters for data
% pre-processing (to get the distance from the target to the measurement
% plane)
% parfile_forprob: a parameter file for inversion
% parfile_vis: a parameter file for visualization. 
% object_properties_func: Matlab function file containing the true properties of the objects
% ObjIdx: a vector of the indices of the targets. Note that all targets
% must have the same distance to the measurement plane in the inversion
% domain. Otherwise, run them separately. 
% thickness_threshold: for Larisa's results only.
% target_depth_file: a matlab file containing the estimated depths of burial. 
% Nguyen Trung Thanh, 2014.
%
% NOTE: the coefficients should be saved in files named "object1.m",
% "object2.m", etc., in the subfolder "final_results" of the main datafolder 
% *The maximum z value in the parameter files for inversion and for
% visualization must be the same. 
% For Larisa test: maybe the step of truncating the coefficient values
% shallower than the burial depth is not needed.

% Warning: do not check the consistency of the input data



% the inversion FEM domain
[x_inv,y_inv,z_inv] = get_FEM_domain(parfile_forprob); 

% the FEM in which we visualize the estimated coeffcient
[x_vis,y_vis,z_vis] = get_FEM_domain(parfile_vis); 

depth_inv = Rx_target_distance_inversion(parfile_preprocessing); % distance from target to measurement plane in the inversion domain


Nobj = length(ObjIdx); % number of objects considered
objprop = object_properties_func; % get the true shapes' properties

dz_vis = z_vis(2) - z_vis(1);

% check the consistency of the parameter files
if abs(z_inv(end) - z_vis(end)) > dz_vis/3
    error('The z-coordinates in the parameter files are not consistent. Correct the parameter file for visualization!');
end


Nx_inv = length(x_inv); 
Ny_inv = length(y_inv); 
Nz_inv = length(z_inv); 

Nx_vis = length(x_vis); 
Ny_vis = length(y_vis); 
Nz_vis = length(z_vis); 


eval(['cd ' datafolder]);

for i = 1:Nobj
    objidx = ObjIdx(i);
    disp(['Object #',num2str(objidx)]);
    CurrentObjProp = objprop(objidx);
    
    if CurrentObjProp.Nrobj==1

        % ------------------------
        % 1. Estimate the size of target from the data: 
        fname = ['object_',num2str(objidx),'/',datafilename]; 
        u = dlmread(fname);  % load the data
        u = reshape(u(1,:),Nx_inv,Ny_inv); % just consider the data at one pseudofrequency 

        
        %interpolate the data from the inversion domain to the visualization domain: 
        [meshX_inv,meshY_inv] = meshgrid(x_inv,y_inv);
        [meshX_vis,meshY_vis] = meshgrid(x_vis,y_vis);

        u_vis = interp2(meshX_inv,meshY_inv,u,meshX_vis,meshY_vis);        
        
        MaskData = estimate_cross_section(abs(u_vis),size_threshold);
        target_size = sum(sum(MaskData)); % size of 
        
        % --------------------------
        % 2. Truncate the reconstructed coefficient:
        fname = [final_result_folder,'/',objfilename,num2str(objidx),'.m'];    
        coef = dlmread(fname); % load the reconstructed coefficient value
        coef = reshape(coef,Nx_inv,Ny_inv,Nz_inv); % convert the coefficient vector to 3-d
        
        %interpolate the coefficient: 
        [meshX_inv,meshY_inv,meshZ_inv] = meshgrid(x_inv,y_inv,z_inv);
        [meshX_vis,meshY_vis,meshZ_vis] = meshgrid(x_vis,y_vis,z_vis);

        coef_vis = interp3(meshX_inv,meshY_inv,meshZ_inv,coef,meshX_vis,meshY_vis,meshZ_vis);     
               
        
        % shift in depth (if needed) and truncate the coefficient before that depth:        
        coef_vis2 = coef_vis*0 + 1; coef_vis = coef_vis(:,:,end:-1:1);                
       
        if nargin > 10
            objdepth = target_depth_file; % load the (true or estimated) burial depths of target:
            if objdepth(objidx,1) == objidx
                burialdepth = objdepth(objidx,2); 
            else
                burialdepth = depth_inv;
            end
        end
        
        shiftIdx = round((burialdepth - depth_inv)/dz_vis);     

        if shiftIdx > 0
            coef_vis2(:,:,1+shiftIdx:end) = coef_vis(:,:,1:end-shiftIdx);
        else
            coef_vis2(:,:,1:end+shiftIdx) = coef_vis(:,:,1-shiftIdx:end);
        end
    
        % truncate the coefficient shallower than the estimated depth: maybe
        % irgnored for Larisa's test::::::::
        depthIdx = round(burialdepth/dz_vis);
        coef_vis2(:,:,1:depthIdx) = 1;
        coef_vis = coef_vis2(:,:,end:-1:1); 
       
        % find the truncation value for size and estimate the cross-section: 
        max_coef = max(max(max(coef_vis)));
        
        if max_coef > 1 % strong target:
            [Trunc_value,Threshold_Coef] = estimate_truncation_value(coef_vis,target_size);
            mask_coef = estimate_cross_section(coef_vis,Threshold_Coef);
        else % weak target:
            [Trunc_value,Threshold_Coef] = estimate_truncation_value(-coef_vis,target_size);
            Trunc_value = -Trunc_value;
            mask_coef = estimate_cross_section(-coef_vis,Threshold_Coef);
        end
        
        % plot the cross section: 
        region_prop = plot_cross_section(MaskData',objidx,x_vis,y_vis,CurrentObjProp); % in visualization, x axis is the column. In data, x axis is the row
        fprintf('%s%5.3f%s%5.3f\n','Centroid in xy plane: ',region_prop.Centroid(1),'  ',region_prop.Centroid(2));
        fprintf('%s%5.3f%s%5.3f\n','Size in x and y directions: ',region_prop.sizex,'  ',region_prop.sizey);
        fprintf('%s%5.3f\n','Truncate threshold for size: ', Trunc_value);
        
        % save the figure of the cross section:
        figname = ['obj',num2str(objidx),'_xy_cs.eps'];
        print('-depsc2',figname);      

        
        % estimate the thickness in the z-direction:
        coef2 = coef_vis; 
        if max_coef > 1
            coef2(coef2 < Trunc_value) = 0;
        else
            coef2(coef2 > Trunc_value) = 0; 
        end
        coef2 = sum(coef2,1); coef2 = sum(coef2,2); coef2(coef2 > 0) = 1;
        disp(['Thickness = ',num2str((sum(coef2,3)-1)*dz_vis)]);

               
        % truncate the coefficient for visualization: 
        % option 1: truncate based on the estimated size:        
%         coef_vis = reshape(coef_vis,1,Nx_vis*Ny_vis*Nz_vis)*bckgr_coef; % without cutoff the coefficient
%         if max_coef > 1
%             coef_vis(coef_vis < Trunc_value) = bckgr_coef;
%             fprintf('%s%5.3f%s%5.3f\n','MaxCoef = ',max_coef,', Threshold for thickness (Larisa test) = ',max_coef*thickness_threshold);
%         else
%             coef_vis(coef_vis > Trunc_value) = bckgr_coef;
%             fprintf('%s%5.3f%s%5.3f\n','MinCoef = ',min(coef_vis),', Threshold for thickness (Larisa test) = ',min(coef_vis)*thickness_threshold);
%         end

        % option 2: truncate based on the cross-section estimated from the
        % data:
        coef_vis = truncate_coef_mask(coef_vis,MaskData);
        coef_vis = reshape(coef_vis,1,Nx_vis*Ny_vis*Nz_vis)*bckgr_coef; % without cutoff the coefficient

        fname = [final_result_folder,'/',objfilename,num2str(objidx),'_trunc.m'];    
        dlmwrite(fname,coef_vis,'delimiter',' ');        
        % Update the parameter files for visualization:             
        ObjCoordinate = get_object_coordinate(region_prop,CurrentObjProp,z_vis(end)-burialdepth);            
        add_target_to_parameter_file(parfile_vis,CurrentObjProp.type,objidx,ObjCoordinate);
    end
    
end

% treate data sets with two targets: 
for i = 1:Nobj
    objidx = ObjIdx(i);
    CurrentObjProp = objprop(objidx);

    if CurrentObjProp.Nrobj==2
        SubObj1 = CurrentObjProp.subobjects(1);
        SubObj2 = CurrentObjProp.subobjects(2);
        
        file1 = [parfile_vis(1:end-4),'_obj',num2str(SubObj1),'.dat'];
        file2 = [parfile_vis(1:end-4),'_obj',num2str(SubObj2),'.dat'];
        
        % modify the parameter file: 
        add_two_targets_to_parameter_file(parfile_vis, objidx, file1, file2);

        % combine the two truncated coefficient file: 
        fname = [final_result_folder,'/',objfilename,num2str(objidx),'_trunc.m'];    
        fname1 = [final_result_folder,'/',objfilename,num2str(SubObj1),'_trunc.m'];    
        fname2 = [final_result_folder,'/',objfilename,num2str(SubObj2),'_trunc.m'];    
        
        scale_for_visualization = 1; % = 0: no scaling, otherwise, we scale for better visualization
        combine_two_coefficient_files(fname1,fname2,fname,bckgr_coef,scale_for_visualization);
        
    end
end
        



close all
