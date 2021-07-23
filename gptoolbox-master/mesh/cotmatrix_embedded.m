function L = cotmatrix_embedded(vers, tris)
  % 返回余切laplace矩阵
 
  % Inputs:
  %   V  #V x dim matrix of vertex coordinates
  %   F  #F x 3  matrix of indices of triangle corners
  % Outputs:
  %   L  sparse #V x #V matrix of cot weights 
 
   % 三角片的边长矩阵。
   edgeLens = [ ...
     sqrt(sum((vers(tris(:,2),:)-vers(tris(:,3),:)).^2,2)) ...      % 第1个顶点的对边
     sqrt(sum((vers(tris(:,3),:)-vers(tris(:,1),:)).^2,2)) ...      % 第2个顶点的对边
     sqrt(sum((vers(tris(:,1),:)-vers(tris(:,2),:)).^2,2)) ...
     ];
 
   L = cotmatrix_intrinsic(edgeLens, tris, size(vers,1));
end

