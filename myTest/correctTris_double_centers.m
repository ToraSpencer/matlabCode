function [ newTris ] = correctTris_double_centers( vers, oldTris, zdir )
newTris = oldTris;
zmax = -inf;
zmin = inf;
versCount = size(vers, 1);
zValue = zeros(versCount, 1);
for i = 1 : versCount
    zValue(i, 1) = dot(vers(i,:), zdir); 
    if zValue(i, 1) > zmax
        zmax = zValue;
    end
    
    if zValue(i, 1) < zmin
        zmin = zValue;
    end
end

middle = (zmax + zmin)/2.0;
versFlag = zeros(versCount, 1);
vers1 = [];
vers2 = [];
for i = 1 : versCount
    if zValue(i, 1) > middle
        versFlag(i, 1) = 1;
        vers2 = [vers2; vers(i, :)];
    else
        vers1 = [vers1; vers(i, :)];
    end
end
center1 = mean(vers1);
center2 = mean(vers2);

trisCount = size(oldTris, 1);
for i = 1: trisCount
    index1 = oldTris(i, 1);
    index2 = oldTris(i, 2);
    index3 = oldTris(i, 3);
    centerToTri = zeros(1, 3);
    if versFlag(index1, 1) == 0 && versFlag(index2, 1) == 0 || ...
            versFlag(index1, 1) == 0 && versFlag(index3, 1) == 0 || ...
            versFlag(index2, 1) == 0 && versFlag(index3, 1) == 0
        centerToTri = vers(index1, :) - center1;
    else
        centerToTri = vers(index1, :) - center2;
    end
    
    norm = cross(vers(index2, :) - vers(index1, :), vers(index3, :) - vers(index1, :));
     
    if dot(centerToTri, norm) <= 0
         temp = newTris(i, 2);
         newTris(i, 2) = newTris(i, 3);
         newTris(i, 3) = temp;
     end
        
end


    

    
end

