function [k,H,K,M,T] = discrete_curvatures(V,F)
%   计算网格的平均曲率，高斯曲率，主曲率；
 
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  % Outputs:
  %   k  #V 主曲率，二元素向量，最大曲率和最小曲率
  %   H  #V list of integrated mean-curvature values
  %   K  #V list of integrated Gaussian curvature

  K = discrete_gaussian_curvature(V,F);
  H = discrete_mean_curvature(V,F);

  M = massmatrix(V,F);

  k = H;

  if nargout > 4
    [A,C] = adjacency_dihedral_angle_matrix(V,F);
    [AI,AJ,AV] = find(A);
    [CI,CJ,CV] = find(C);
    assert(isequal(CI,AI));
    assert(isequal(CJ,AJ));

    
    opp = sub2ind(size(l),CI,CV);
    inc = [ ...
      sub2ind(size(l),CI,mod(CV+1-1,3)+1) ...
      sub2ind(size(l),CI,mod(CV+2-1,3)+1)];
    lV = l(opp);

    T = zeros(size(V,1),3,2);
    for v = 1:size(V,1)
    end
  end

end
