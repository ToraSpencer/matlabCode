function [ Bcenter ] = barycenter(vers, tris)
 
  % 计算每个三角片的重心
  % Inputs:
  %   V #V x 3 matrix of vertex coordinates
  %   F #F x 3  matrix of indices of triangle corners
  % Output:
  %   B a #F x 3 matrix of 3d vertices

  i1 = tris(:,1);
  i2 = tris(:,2);
  i3 = tris(:,3);
  
  Bcenter = (vers(i1,:) + vers(i2,:) + vers(i3,:)) / 3.0;

end

