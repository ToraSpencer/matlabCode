function pts = critical_points(V,F,u,tolerance)

%REVERSE_SUBDIVISION 
% Find out which face contains the critical points of a
 
% Inputs:
%  V,F  the coarse input mesh
%  u  the function for which to find critical points
%  tolerance  the tolerance for writical points (if the gradient norm is below this,
%                                          we are at a critical point)
% Outputs:
%  pts  a list of face indices into the rows of F that tell us which faces
%       contain the critical points.
%

G = grad(V,F);
gu = G*u;
pts = find(normrow(gu) < tolerance);

end

