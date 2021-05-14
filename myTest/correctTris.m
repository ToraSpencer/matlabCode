function [ newTris ] = correctTris(vers, oldTris)
 center = mean(vers);
 triCount = size(oldTris, 1);
 newTris = oldTris;
 
 for i = 1:triCount
     currentTri = oldTris(i, :);
     index1 = currentTri(1);
     index2 = currentTri(2);
     index3 = currentTri(3);
     centerToTri = vers(index1, :) - center;
     norm = cross(vers(index2, :) - vers(index1, :), vers(index3, :) - vers(index1, :));
     if dot(centerToTri, norm) <= 0
         temp = newTris(i, 2);
         newTris(i, 2) = newTris(i, 3);
         newTris(i, 3) = temp;
     end
 end


end

