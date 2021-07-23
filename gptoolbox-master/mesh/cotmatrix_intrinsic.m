function L = cotmatrix_intrinsic(edgeLens, tris, versCount)
  % 使用三角片及其边长数据计算余切laplacian
 
  % Inputs:
  %   edgeLens  #F by 3, array of edge lengths of edges opposite each face in F
  %   F  #F by 3, list of indices of triangle corners
  %   versCount  number of vertices, only needed to set size
  % Outputs:
  %   L  sparse nvert x nvert matrix of cot weights 

  if(size(tris,1) == 3)
    warning('F seems to be 3 by #F, it should be #F by 3');
  end
  tris = tris';

  % Law of cosines + Law of sine gives you:

  % renaming indices of vertices of triangles for convenience
  i1 = tris(1,:); i2 = tris(2,:); i3 = tris(3,:); 
  l1 = edgeLens(:,1); l2 = edgeLens(:,2); l3 = edgeLens(:,3);
  
  % 三角片周长的一半
  s = (l1 + l2 + l3)*0.5;
  
  
  % 使用Heron公式计算三角片面积
  dblA = 2*sqrt( s.*(s-l1).*(s-l2).*(s-l3));          % 三角片面积
  
  % 计算余切值：
  % 余弦公式：cosC = (a^2+b^2-c^2)/(2*a*b);
  % 三角片面积：dblA = 1/2*a*b*sinC   →   sinC = 2*dblA/(a*b);
  % 余切：cotC = cosC/sinC = (a^2+b^2-c^2)/(4*dblA);
  cot3 = (l1.^2 + l2.^2 -l3.^2)./dblA/4; 
  cot1 = (l2.^2 + l3.^2 -l1.^2)./dblA/4; 
  cot2 = (l1.^2 + l3.^2 -l2.^2)./dblA/4; 
  
  % diag entries computed from the condition that rows of the matrix sum up to 1
  % (follows from  the element matrix formula E_{ij} = (v_i dot v_j)/4/A )
  diag1 = -cot3-cot2; diag2 = -cot3-cot1; diag3 = -cot2-cot1;
  % indices of nonzero elements in the matrix for sparse() constructor
  
  i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
  j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
  
  % values corresponding to pairs form (i,j)
  v = [cot3 cot3 cot1 cot1 cot2 cot2 diag1 diag2 diag3];
  
  % for repeated indices (i,j) sparse automatically sums up elements, as we want
  L = sparse(i, j, v, versCount, versCount);
end
