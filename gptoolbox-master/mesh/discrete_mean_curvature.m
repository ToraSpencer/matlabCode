function [H] = discrete_mean_curvature(vers, tris)
  % 计算网格每个顶点的平均曲率；

  
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  % Outputs:
  %   H  #V list of integrated mean-curvature values

  
  % Examples:
  %   H = discrete_mean_curvature(V,F);
  %   M = massmatrix(V,F);
  %   tsurf(F,V,'CData',M\H,'EdgeColor','none',fphong);
 
  [adjAngles, C] = adjacency_dihedral_angle_matrix(vers,tris);
  [AI,AJ,AV] = find(adjAngles);
  [CI,CJ,CV] = find(C);
  assert(isequal(CI,AI));
  assert(isequal(CJ,AJ));
  l = edge_lengths(vers,tris);
 
  opp = sub2ind(size(l),CI,CV);
  inc = [ ...
    sub2ind(size(l),CI,mod(CV+1-1,3)+1) ...
    sub2ind(size(l),CI,mod(CV+2-1,3)+1)];
  lV = l(opp);
 
  H = full(sparse(tris(inc),1,0.5*0.5*0.5*repmat((pi-AV).*l(opp),1,2),size(vers,1),1));
end
