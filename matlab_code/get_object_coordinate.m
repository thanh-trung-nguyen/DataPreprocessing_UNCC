function ObjCoord  = get_object_coordinate(region_prop,objprop,burialdepth)
% get the coordinates of an object from the regional properties and true
% object properties. This is to be added to the parameter file for
% visualization.

ObjCoord = zeros(1, 6);

if objprop.Nrobj == 1    
    ObjType = objprop.type;
    if strcmpi(ObjType,'prism')

        x1 = region_prop.Centroid(1) - objprop.sizex/2;
        x2 = region_prop.Centroid(1) + objprop.sizex/2;
        y1 = region_prop.Centroid(2) - objprop.sizey/2;
        y2 = region_prop.Centroid(2) + objprop.sizey/2;

        ObjCoord = [x1, x2, y1, y2, burialdepth-objprop.sizez, burialdepth];

    elseif strcmpi(ObjType,'sphere') 
        centerX = region_prop.Centroid(1);
        centerY = region_prop.Centroid(2);
        centerZ = burialdepth - objprop.radius;
        ObjCoord = [centerX, centerY, centerZ, objprop.radius, -9999 -9999];
        
    elseif strcmpi(ObjType,'cylinder-z')
        centerX = region_prop.Centroid(1);
        centerY = region_prop.Centroid(2);
        
        ObjCoord = [centerX,c centerY, objprop.radius, burialdepth-objprop.sizez, burialdepth];
                
    elseif strcmpi(ObjType,'cylinder-x')
        x1 = region_prop.Centroid(1) - objprop.sizex/2;
        x2 = region_prop.Centroid(1) + objprop.sizex/2;
   
        centerY = region_prop.Centroid(2);
        centerZ = burialdepth - objprop.radius;
        
        ObjCoord = [centerY, centerZ, objprop.radius, x1, x2, -9999];
        
    elseif strcmpi(ObjType,'cylinder-y')

        centerX = region_prop.Centroid(1);
        centerZ = burialdepth - objprop.radius;
        y1 = region_prop.Centroid(2) - objprop.sizey/2;
        y2 = region_prop.Centroid(2) + objprop.sizey/2;

        ObjCoord = [centerX, centerZ, objprop.radius, y1, y2, -9999];
        
        
    elseif strcmpi(ObjType,'matreoshka');
        disp('Matreoshka: manual updating the parameter file needed');
        
    end

elseif objprop.Nrobj == 2    
    disp('Two targets: copy the object properties from the two halves!');

else
    disp('Only 1 or 2 objects allowed');
end
