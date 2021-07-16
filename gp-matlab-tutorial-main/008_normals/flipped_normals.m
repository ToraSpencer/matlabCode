function FN = flipped_normals(vers, tris)
% Compute the flipped normals of a triangle mesh
% flipped normal——翻转的法线方向，即FN = -normalDir;
 

%Compute per-face normals.
N = normals(vers,tris);
N = N ./ normrow(vers,tris);

%Flip the per-face normals.
FN = -N;

end

