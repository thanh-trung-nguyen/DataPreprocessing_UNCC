function stat = plot_cross_section(MaskFile,objidx,x,y,objprop)
% display the cross section and the ground true (if provided)
% the output "stat" is the properties of the region of th target in the xy
% coordinate.
% Nguyen Trung Thanh, 2014.
% NOTE that x is the horizontal coordinate, y is the vertical coordinate

[Ny,Nx] = size(MaskFile);

% Plot the cross-section and the ground truth: 
u3 = MaskFile;
u3(u3 > 0) = 30; u3(u3==0) = 63;
figure; set(gca,'fontsize',20);
image(u3); colormap(gray); truesize(size(u3)*8); 

if nargin > 4
    
    Tick_x = 1:20:Nx;
    Tick_y = 1:20:Ny;
    dx = x(2) - x(1); 
    dy = y(2) - y(1);
    
    MaskFile2 = logical(MaskFile);
    stat = regionprops(MaskFile2); % get properties of regions of the mask file
    % the true object:
    if objprop.Nrobj == 1    
        ObjType = objprop.type;
        if strcmpi(ObjType,'prism')
            x1 = stat.Centroid(1) - objprop.sizex/dx/2;
            x2 = stat.Centroid(1) + objprop.sizex/dx/2;
            y1 = stat.Centroid(2) - objprop.sizey/dy/2;
            y2 = stat.Centroid(2) + objprop.sizey/dy/2;

            X = [x1 x2 x2 x1 x1];
            Y = [y1 y1 y2 y2 y1];
            
            
        elseif strcmpi(ObjType,'sphere') 
            t = 0:0.001:2*pi; 
            X = objprop.radius/dx*cos(t) + stat.Centroid(1);
            Y = objprop.radius/dy*sin(t) + stat.Centroid(2);  
            
        elseif strcmpi(ObjType,'cylinder-z')
            t = 0:0.001:2*pi; 
            X = objprop.radius/dx*cos(t) + stat.Centroid(1);
            Y = objprop.radius/dy*sin(t) + stat.Centroid(2);       

        elseif strcmpi(ObjType,'cylinder-x')
            x1 = stat.Centroid(1) - objprop.sizex/dx/2;
            x2 = stat.Centroid(1) + objprop.sizex/dx/2;
            y1 = stat.Centroid(2) - objprop.radius/dy;
            y2 = stat.Centroid(2) + objprop.radius/dy;

            X = [x1 x2 x2 x1 x1];
            Y = [y1 y1 y2 y2 y1];    

        elseif strcmpi(ObjType,'cylinder-y')
            x1 = stat.Centroid(1) - objprop.radius/dx;
            x2 = stat.Centroid(1) + objprop.radius/dx;
            y1 = stat.Centroid(2) - objprop.sizey/dy/2;
            y2 = stat.Centroid(2) + objprop.sizey/dy/2;

            X = [x1 x2 x2 x1 x1];
            Y = [y1 y1 y2 y2 y1];  

        elseif strcmpi(ObjType,'matreoshka');
            [X,Y] = matreoshka;        
            X = -X/dx + stat.Centroid(1);
            Y = -Y/dy + stat.Centroid(2)+1;
        end
        
        % plot the true cross-section:
        hold on; plot(X,Y,'-k','linewidth',2); hold off;
        set(gca,'xtick',Tick_x,'xticklabel',x(Tick_x));
        set(gca,'ytick',Tick_y,'yticklabel',y(Tick_y));    
        xlabel('x (m)'); ylabel('y (m)');

%         % save for two targets cases: 
%         fname = ['obj',num2str(objidx),'_cross_section.mat'];  
%         eval(['save ' fname ' MaskFile X Y x y Tick_x Tick_y']);

    elseif objprop.Nrobj == 2    
        disp('Two targets simulteneously are under development');

    else
        disp('Only 1 or 2 objects allowed');
    end
end

% convert from image coordinate to the original coordinate:
stat.Centroid(1) = (stat.Centroid(1)-1)*dx + x(1);
stat.Centroid(2) = (stat.Centroid(2)-1)*dy + y(1);
stat.sizex = stat.BoundingBox(3)*dx;
stat.sizey = stat.BoundingBox(4)*dx;


