function [bdryEdges] = outline(tries)
  % 找出三维网格的所有边缘边
 

  %%
  %% This does not maintain original order
  %%
  %% Find all edges in mesh, note internal edges are repeated
  %E = sort([F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)]')';
  %% determine uniqueness of edges
  %[u,m,n] = unique(E,'rows');
  %% determine counts for each unique edge
  %counts = accumarray(n(:), 1);
  
  %% extract edges that only occurred once
  %O = u(counts==1,:);

  % build directed adjacency matrix
  A = sparse(tries,tries(:,[2:end 1]),1);
  
  % Find single occurance edges
  [OI,OJ,OV] = find(A-A');
  
  % Maintain direction
  bdryEdges = [OI(OV>0) OJ(OV>0)];%;OJ(OV<0) OI(OV<0)];

end
