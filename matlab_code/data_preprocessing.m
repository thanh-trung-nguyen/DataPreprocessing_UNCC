function data_preprocessing(datafolder,paramfile,startingstep,endstep,incwavefile,scalingfactorfile,Shift)
% 
% Function "data_preprocessing.m" performs the data preprocessing for UNCC
% data.
% 
% Steps of data preprocessing: 
% Step 1: zero-offset shift
% Step 2: upsampling and resizing data
% Step 3: time-zero correction
% Step 4: propagate the data from the measurement plane to a closer plane
% Step 5: Shift the data to a given location of source as in the inversion domain
% Step 6: extract the target's signals
% Step 7: Scale the data and prepare data for inversion: change domain, compute the Laplace
% transform, ...
%
% USE: 
% To run the file, use the 
% input: datafolder: folder containing the input data file, result of the function
% "loaddata_UNCC.C"
% paramfile: a structured text file containing the necessary parameters for
% preprocessing
% startingstep: step of preprocessing that is started with
% endstep: the last step of data preprocessing you want to run
% invwavefile: this is a simulated solution in homogeneous medium
% scalingfactorfile: a file which contains a vector of scaling factor in
% pseudo-frequency domain.
% Shift: This parameter is provided in case that the first direct signal is
% not detected. Instead, the first signal which can be detected has some
% shift from the true first direct signal. Shift = positive if the detected
% signal is after the true first direct signal. In UNCC data, this number
% is 164. 

% Output: a series of files of the results after certain number of
% preprocessing steps.
%  
% Nguyen Trung Thanh, UNCC 2013
% -------------------------------------------------------------------------

lightspeed = 0.3; % light speed in meter/nanoseconds

if nargin < 7
    Shift = 0;
end

if (nargin < 3)
    startingstep = 1;
    endstep = 8;    
end
if (startingstep > endstep)
    endstep = startingstep;
end

% ----------------file names (new update on Feb 22, 2013: change the file names): 
if strcmp(datafolder(end),'/') == 0
    datafolder = [datafolder,'/'];
end

fname_t0 = [datafolder,'meadat_t0.dat'];
fname_propagated = [datafolder,'prop_dat.dat'];
fname_target = [datafolder,'prop_dat_target.dat'];
fname_source_shift = [datafolder,'prop_dat_source_shift.dat'];
fname_target_scale = [datafolder,'prop_dat_target_scale.dat'];
fname_inversion = [datafolder,'prop_dat_inv.dat'];
fname_inversion_scale = [datafolder,'prop_dat_inv_scale.dat'];

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

Nx_mea = round((Xmax_mea - Xmin_mea)/dx_mea) + 1;
Ny_mea = round((Ymax_mea - Ymin_mea)/dy_mea) + 1;

X_mea = (linspace(Xmin_mea,Xmax_mea,Nx_mea))'; 
Y_mea = (linspace(Ymin_mea,Ymax_mea,Ny_mea))';
Np_mea = Nx_mea*Ny_mea;

Nx_FDM = round((Xmax_FDM - Xmin_FDM)/dx_FDM) + 1;
Ny_FDM = round((Ymax_FDM - Ymin_FDM)/dy_FDM) + 1;

X_FDM = (linspace(Xmin_FDM,Xmax_FDM,Nx_FDM))'; 
Y_FDM = (linspace(Ymin_FDM,Ymax_FDM,Ny_FDM))';
Np_FDM = Nx_FDM*Ny_FDM;

Nx_FEM = round((Xmax_FEM - Xmin_FEM)/dx_FEM) + 1;
Ny_FEM = round((Ymax_FEM - Ymin_FEM)/dy_FEM) + 1;

X_FEM = (linspace(Xmin_FEM,Xmax_FEM,Nx_FEM))'; 
Y_FEM = (linspace(Ymin_FEM,Ymax_FEM,Ny_FEM))';
Np_FEM = Nx_FEM*Ny_FEM;

MaxTime = MaxTime*lightspeed;
dt = dt*lightspeed; % dimensionless time step
dtnew = dt/UpsamplingFactor; 
Nt = round(MaxTime/dt) + 1; % number of data points in time after upsampling
Ntnew = round(MaxTime/dtnew) + 1;

Ns = round((s_max-s_min)/ds) + 1;
s = linspace(s_max,s_min,Ns); % frequencies of Laplace transform


% ========= STEPS OF PREPROCESSING: 

if (startingstep <= 3) && (endstep >= 3)
    disp('Step 1: Offset correction: '); 
    % ------- Step 1: offset correction
    meadatfile = 'measured_data.dat'; 
    if ~exist([datafolder,meadatfile],'file')
        meadatfile = 'measured_data.m'; % for old files
    end
    data = dlmread([datafolder,meadatfile]);
    data = offset_correction(data);

    % ------- Step 2: resizing data and upsampling:    
    if Nt ~= size(data,1);
        disp('Step 2: Resize data in time and upsampling:'); 
        dataNew = zeros(Nt,Np_mea);
        dataNew(1:size(data,1),:) = data; % add zero to the end if needed
    end
    t = linspace(0,MaxTime,Nt);
    data = upsampling(dataNew,t,UpsamplingFactor,'linear');

    % ------- Step 3: time zero estimation:
    disp('Step 3: Time zero correction: ');
    data = time_zero_correction(data,distTxRx,dtnew,NoiseLevel,Shift);
    dlmwrite(fname_t0,data,'delimiter',' '); % write the data to a file
end


% ------ Step 4: Propagate the data:
if (startingstep <= 5) && (endstep >= 4)
    disp('Step 4: Data propagation:'); 
    if ~exist('data','var');
        data = dlmread(fname_t0);
    end
    
    Ncenter = [round((Nx_mea+1)/2), round((Ny_mea+1)/2)];
    minimum_distance = 2*EstDistRxObj - distTxRx - 0.2; 
    dmea = Rx_target_distance(data(:,IdxProp*Nx_mea+Ncenter(1)),minimum_distance,dtnew,1.5*NoiseLevel);
    dist1 = abs(IdxProp - Ncenter(2))*dy_mea;
    distRxObj = find_correct_distance(distTxRx,dmea,dist1) + distTxRx
    PropagationDistEst = distRxObj - distRxObj_inv; % in case the target is in the air, we propagate up to 4 cm distance
    if PropagationDist < 0
        PropagationDist = PropagationDistEst;
    end
    
    tnew = linspace(0,MaxTime,Ntnew);
    data = reshape(data,Ntnew,Nx_mea,Ny_mea);
    N = round(minimum_distance/dtnew); data(1:N,:,:) = 0; % remove cross-talk and unwanted signals earlier than the target's signal
    data = wave_propagation_3d(data,PropagationDist,tnew,X_mea,Y_mea);
    data = reshape(data,Ntnew,Np_mea); 
    fprintf('%s%d\n','Minimum value of the propagated data: ',min(min(data(1:end,:))));
    dlmwrite(fname_propagated,data,'delimiter',' ');

    disp('Step 5: Source shift:'); 

    distTxObj = distRxObj -distTxRx;
    data = source_shift(data,distTxObj_inv,distRxObj_inv,dtnew,distTxObj);    
%     data = source_shift(data,distTxObj_inv,distRxObj_inv,dtnew);    
    dlmwrite(fname_source_shift,data,'delimiter',' ');    

end



%%%%% July 15, 2014: NEED UPDATE from this step.
% ------ Step 6: extract the target's signals: 
if (startingstep <= 6) && (endstep >=6)
    disp('Step 6: Extract target signals:'); 
    if ~exist('data','var');
        data = dlmread(fname_source_shift);
    end
    data = extract_target_signal(data);
    dlmwrite(fname_target,data,'delimiter',' ');    
end
    
% ------ Step 7: scaling and prepare data for inversion: 
if (startingstep <= 7) && (nargin > 4) && (endstep >=7) % ignored if the incident wave is not provided
    % ------ Step 7: scaling: 
    if nargin > 5
        scalingfactor = dlmread(scalingfactorfile); % load the scaling factor which has been estimated in advance for Laplace transform
%         scalingfactor = scalingfactor*1.1; % compensate for the effect of sand
    else
        scalingfactor = 0*s + 1;
    end

    ScalingFactor = mean(scalingfactor);  % scaling factor for the time domain data.

    
    disp('Step 7: Prepare data for the GCA: '); 
    if ~exist('data','var');
        data = dlmread(fname_target);
    end
    
    % resizing time-domain data to the FEM inversion domain:
    disp('Resizing time domain data to the FEM inversion domain:');
    dataNew = zeros(Nt_inv,Np_FEM);
    for idt = 2:Nt_inv+1
        u = reshape(data(idt,:),Nx_mea,Ny_mea); 
        u = bilinear_interpolation(u,X_mea,Y_mea,X_FEM,Y_FEM);
        dataNew(idt-1,:) = reshape(u,1,Np_FEM);
    end
    data = dataNew;
    dlmwrite(fname_inversion,data,'delimiter',' '); % time domain data before scaling
    dlmwrite(fname_inversion_scale,data*ScalingFactor,'delimiter',' '); % time domain data after scaling
    
    %find the strongest signal:
    MinPeak = min(min(data));
    [n,MinIdx] = find(data==MinPeak);
    MinIdx = MinIdx(1);
    
    % compute the Laplace transform of the data and average:
    disp('Compute Laplace transform of the data:');
    ws = laplacetransform(data,dtnew,s);
    ws(ws > 0) = 0; % remove "unusual" data points.
    
    
    % averaging and scaling: 
    for k = 1:Ns
          w2 = ws(k,:);                   
          MinValue = ws(k,MinIdx);
     
          Idx1 = find(w2 < MinValue);
          w2(Idx1) = MinValue*min(data(:,Idx1))/MinPeak;% replace values smaller than that of the strongest signal

          for m = 1:5 % averaging but keep the minimum value
                w2 = average(w2,Nx_FEM,Ny_FEM); w2 = w2*MinValue/min(min(w2));
          end
        ws(k,:) = w2*scalingfactor(k); 
    end 
        
    % compute Laplace transform of the incident wave:
    disp('Compute Laplace transform of incident wave if not yet done:');
    idx = find(incwavefile == '.'); 
    fname1 = sprintf('%s%s%3.2f%s%3.2f%s%d%s',incwavefile(1:idx-1),'_w_s',s(1),'_',s(Ns),'_',Ns,incwavefile(idx:end));    
    file1 = dir(incwavefile); date1 = file1.datenum; 

    ReadFromFile = 0; 
    if exist(fname1,'file')
        file2 = dir(fname1); date2 = file2.datenum;     
        if date1 < date2
            ReadFromFile = 1;
        end
    end        
    if ~ReadFromFile % if not present or out of date, recompute
        Ui = dlmread(incwavefile);
        Wi = laplacetransform(Ui,dtnew,s);
        dlmwrite(fname1,Wi,'delimiter',' ');
    else
        Wi = dlmread(fname1); % load if available, so only compute once!
    end
    % interpolate the incident wave to the FEM inversion domain:
    wi = zeros(Ns,Np_FEM);
    for k = 1:Ns
        W0 = reshape(Wi(k,:),Nx_FDM,Ny_FDM);
        wi(k,:) = reshape(bilinear_interpolation(W0,X_FDM,Y_FDM,X_FEM,Y_FEM),1,Np_FEM);
    end
        
    % compute function v = log(1 + ws/wi)/s^2:
    disp('Compute function v and save the result: ');
    v_gca = zeros(Ns,Np_FEM);
    for k = 1:Ns
        v_gca(k,:) = log(1 + ws(k,:)./wi(k,:))/s(k)^2;  
    end
    fname_v = sprintf('%s%s%3.2f%s%3.2f%s%d%s',datafolder,'prop_dat_inv_scale_v_s',s(1),'_',s(Ns),'_',Ns,'.m');    
    dlmwrite(fname_v,v_gca,'delimiter',' ');
    
    
    
    % ---------------Create input file for Thanh code:
    % get the boundary grid information: 
    idx = find(incwavefile == '/',1,'last'); 
    boundarygridfile = [incwavefile(1:idx),'grid_boundary.inp'];
    fid = fopen(boundarygridfile);
    if fid == -1
        error('The boundary grid file not found');
    end
    NoBndNodes = round(str2double(fscanf(fid,'%s',1))); fprintf('%s%d\n','Number of boundary nodes: ',NoBndNodes); 

    data = zeros(Ns+1,NoBndNodes); 
    
    xb = zeros(1,NoBndNodes); yb = xb; zb = xb; 
    for i = 1:NoBndNodes
        data(1,i) = round(str2double(fscanf(fid,'%s',1)));
        xb(i) = str2double(fscanf(fid,'%s',1));
        yb(i) = str2double(fscanf(fid,'%s',1));
        zb(i) = str2double(fscanf(fid,'%s',1));
    end
    fclose(fid);    
%     % find locations of measured data base on the coordinate fit: TEST 1
%     for nb = 1:NoBndNodes
%         if (abs(zb(nb) - z_inv) < 1e-6)
%             for idx = 1:Np1
%                if (abs(X1(idx) - xb(nb))< 1e-6 && abs(Y1(idx) - yb(nb)) < 1e-6)
%                     data(2:end,nb) = v_gca(:,idx);
%                end
%             end
%         end
%     end
    
    data(2:end,end-Np_FEM+1:end) = v_gca; % data at the backscattering side TEST 2
  

    % computing q: 
    fname_v0 = sprintf('%s%s%3.2f%s%3.2f%s%d%s',incwavefile(1:idx),'SolFEMBoundary_v_s',s(1),'_',s(Ns),'_',Ns,'.m');    
    fname_solbnd = [incwavefile(1:idx),'SolFEMBoundary.m'];
    file1 = dir(fname_solbnd); date1 = file1.datenum; 

    ReadFromFile = 0;
    if exist(fname_v0,'file')
        file2 = dir(fname_v0); date2 = file2.datenum;     
        if (date1 < date2)
            ReadFromFile = 1;
        end
    end
    if ~ReadFromFile
        u0 = dlmread(fname_solbnd);
        v0 = compute_v(u0(2:end,:),dtnew,s); 
        dlmwrite(fname_v0,v0,'delimiter',' ');
    else
        v0 = dlmread(fname_v0);
    end
    
    v = data(2:end,:) + v0;
    data(2:Ns,:) = (v(1:Ns-1,:) - v(2:Ns,:))/(s(1)-s(2)); % function q
    data(Ns+1,:) = v(1,:); % the tail
   
    bnddatafile= sprintf('%s%s%3.2f%s%3.2f%s%d%s',datafolder,'boundarydata_s',s(1),'_',s(Ns),'_',Ns-1,'.m');
    save_file_mix_data(bnddatafile,data(1,:), data(2:end,:));
    
    
    
    % =====================================================================
    %-------create the masks for stopping criterion (for Thanh's code):
    V = abs(reshape(v_gca(1,:),Nx_FEM,Ny_FEM));    
    
%     V(2:end-1,2:end-1) = V(1:end-2,2:end-1) + V(2:end-1,2:end-1) + V(3:end,2:end-1) + V(2:end-1,1:end-2) + V(2:end-1,3:end);
    V(V < 0.65*max(max(V))) = 0;     V(V> 0) = 1;
      
    V2 = V;    V2(V2 > 0) = 1;

    dlmwrite([datafolder,'/coef_mask.m'],V2','delimiter',' ');

    figure; imagesc(V2'); truesize([Ny_FEM,Nx_FEM]*10);
    
    figname = [datafolder,'/coef_mask.jpg'];      
    print('-djpeg90', figname);

    w3 = mean(V,2); 
    w4 = mean(V,1);

    idx1 = find(w3 > 0);
    idx2 = find(w4 > 0);
    NrIdx = 1;
    V2(idx1(1)-NrIdx:idx1(end)+NrIdx,idx2(1)-NrIdx:idx2(end)+NrIdx) = 1;

    dlmwrite([datafolder,'/coef_mask_large.m'],V2','delimiter',' ');

    figure; imagesc(V2'); truesize([Ny_FEM,Nx_FEM]*10);
    figname = [datafolder,'/coef_mask_large.jpg'];      
    print('-djpeg90', figname);       


end

