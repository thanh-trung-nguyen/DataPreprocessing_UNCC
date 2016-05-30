function depth = Rx_target_distance_inversion(paramfile)
% get the distance from the measurement plane to the target in the
% inversion domain. This parameter is given in the parameter file for data
% preprocessing. 


% -----------------load the parameters: 
fid = fopen(paramfile);
if fid == -1
    error('The parameter file not found');
end
text = fscanf(fid,'%s',1);  Target_is_Buried = fscanf(fid,'%s',1); 
text = fscanf(fid,'%s',1);  MaxTime = str2double(fscanf(fid,'%s',1)); 
text = fscanf(fid,'%s',1); dt = str2double(fscanf(fid,'%s',1)); %#ok<*NASGU>
text = fscanf(fid,'%s',1); UpsamplingFactor = str2double(fscanf(fid,'%s',1));
text = fscanf(fid,'%s',1); Xmin_mea = str2double(fscanf(fid,'%s',1)); Xmax_mea = str2double(fscanf(fid,'%s',1)); dx_mea = str2double(fscanf(fid,'%s',1)); 
text = fscanf(fid,'%s',1); Ymin_mea = str2double(fscanf(fid,'%s',1)); Ymax_mea = str2double(fscanf(fid,'%s',1)); dy_mea = str2double(fscanf(fid,'%s',1)); 
text = fscanf(fid,'%s',1); distTxRx = str2double(fscanf(fid,'%s',1));
text = fscanf(fid,'%s',1); EstDistRxObj = str2double(fscanf(fid,'%s',1));
text = fscanf(fid,'%s',1); NoiseLevel = str2double(fscanf(fid,'%s',1));
text = fscanf(fid,'%s',1); ScalingFactor = str2double(fscanf(fid,'%s',1));
text = fscanf(fid,'%s',1); PropagationDist = str2double(fscanf(fid,'%s',1)); IdxProp = round(str2double(fscanf(fid,'%s',1)));
text = fscanf(fid,'%s',1); Nt_inv = round(str2double(fscanf(fid,'%s',1)));
text = fscanf(fid,'%s',1); distTxObj_inv = str2double(fscanf(fid,'%s',1));
text = fscanf(fid,'%s',1); distRxObj_inv = str2double(fscanf(fid,'%s',1));
text = fscanf(fid,'%s',1); Xmin_FDM = str2double(fscanf(fid,'%s',1)); Xmax_FDM = str2double(fscanf(fid,'%s',1)); dx_FDM = str2double(fscanf(fid,'%s',1)); 
text = fscanf(fid,'%s',1); Ymin_FDM = str2double(fscanf(fid,'%s',1)); Ymax_FDM = str2double(fscanf(fid,'%s',1)); dy_FDM = str2double(fscanf(fid,'%s',1)); 
text = fscanf(fid,'%s',1); Xmin_FEM = str2double(fscanf(fid,'%s',1)); Xmax_FEM = str2double(fscanf(fid,'%s',1)); dx_FEM = str2double(fscanf(fid,'%s',1)); 
text = fscanf(fid,'%s',1); Ymin_FEM = str2double(fscanf(fid,'%s',1)); Ymax_FEM = str2double(fscanf(fid,'%s',1)); dy_FEM = str2double(fscanf(fid,'%s',1)); 
text = fscanf(fid,'%s',1); s_min = str2double(fscanf(fid,'%s',1)); s_max = str2double(fscanf(fid,'%s',1)); ds = str2double(fscanf(fid,'%s',1)); 
fclose(fid);

depth = distRxObj_inv;
