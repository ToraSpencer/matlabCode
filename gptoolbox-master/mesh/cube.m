function [vers, tris, Q] = cube(x,y,z)
  % CUBE Construct a mesh of the unit cube. 
  % Sides are ordered like sides of a die (one of many dice).
  % Inputs:
  %   x  number of vertices along x-axis
  %   y  number of vertices along y-ayis
  %   z  number of vertices along z-azis
  % Outputs:
  %   V  x*y*z by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  %   Q  #Q by 3 list of quad indices

  if nargin<2
    y = x;
  end
  if nargin<3
    z = y;
  end

  sam = [x y;z y;x z;x z;z y;x y];
  axes = [0 1 0; 0 1 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0];
  angles = [0 pi/2 pi/2 -pi/2 -pi/2 pi];
  vers = [];
  tris = [];
  
  for s = 1:6
    [CV,CF] = create_regular_grid(sam(s,1),sam(s,2),0,0);
    CV(:,3) = 0;
    R = round(axisangle2matrix(axes(s,:),angles(s)));
    tris = [tris;size(vers,1)+CF];
    vers = [vers;(CV-0.5)*R+0.5];
  end
  
  Q = [tris(1:2:end-1,[1 2]) tris(2:2:end,[2 3])];
  
  % Should be able to do this procedurally
  [vers,~,J] = remove_duplicate_vertices(vers,1e-12);
  Q = J(Q);
  tris = J(tris);
  
  % oops inside out
  tris = fliplr(tris);              % ·­×ªÈý½ÇÆ¬£»
  Q = fliplr(Q);

end
