function combine_two_coefficient_files(inpfile1,inpfile2,outfile,bckgr_coef,scale_for_visualization)
% add results of two halves data into the complete data sets. 
% example: 
% inpfile1 = 'obj16_trunc.m';
% inpfile2 = 'obj17_trunc.m';
% outfile = 'obj9_trunc.m';
% scale_for_visualization: 1 if yes, 0 if No. defaut is No. 

if nargin < 5
    scale_for_visualization = 0;
end

u1 = dlmread(inpfile1);
u2 = dlmread(inpfile2);

max1 = max(u1); min1 = min(u1); 
max2 = max(u2); min2 = min(u2);

if max1 > bckgr_coef && max2 > bckgr_coef % strong targets:	  
	u1(u1 <= bckgr_coef) = 0;
    if scale_for_visualization ~= 0
        u1 = u1*max2/max1; % scale the first coefficient to the same level se the second one, for better visualization in AVS viewer
    end
	u = max(u1,u2);
elseif max1 < bckgr_coef && max2 < bck_coef % weak targets:
	u1(u1 >= bckgr_coef) = 0; 
    if scale_for_visualization ~=0
        u1 = u1*min2/min1; % scale the first coefficient to the same level se the second one, for better visualization in AVS viewer
    end
	u = u1 + u2;
end

dlmwrite(outfile,u,'delimiter',' ');
