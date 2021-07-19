function [adjAngles, C] = adjacency_dihedral_angle_matrix(vers, tris)
% 生成相邻二面角矩阵
% 若三角片i和三角片j有公共边，则A(i, j) == 相邻二面角弧度；若没有公共边则元素为0；

  % Outputs:
  %   adjAngles  #F by #F sparse matrix of signed dihedral angles. All entries are in
  %     [0,2*pi] adjAngles-->0 as edge is more convex, adjAngles-->pi as edge is more flat,
  %     adjAngles-->2*pi as edge is more concave (as view from outside, assuming
  %     counter-clockwise orientation.
  %     二面角是锐角还是钝角是从外部看起来
  
  %   C  #F by #F sparse matrix revealing C(i,j) = c that face j is incident on
  %     face i across its corner c

  
  
  % Example:
  %   adjAngles = adjacency_dihedral_angle_matrix(V,F);
  %   % unsigned dihedral angles
  %   UA = pi*(A~=0)-abs(pi*(A~=0)-A);
  % 
  %   adjAngles = adjacency_dihedral_angle_matrix(V,F);
  %   % Adjacency matrix of nearly coplanar neighbors
  %   UA = pi*(A~=0)-abs(pi*(A~=0)-A);
  %   AF = UA>=(pi-1e-5);
  %   % get connected components
  %   C = components(AF);
  %   % plot with unique colors for each coplanar patch
  %   set(tsurf(F,V),'CData',C);
  %   CM = jet(numel(unique(C)));
  %   colormap(CM(randperm(end),:));
  %

  % all edges "both" directions
  edges = [tris(:,[2 3]);tris(:,[3 1]);tris(:,[2 1])];
  
  % index for each edge to unique edge
  [sortedEdges, ~, IC] = unique(sort(edges, 2), 'rows');
  
  % FE(i,j) = 1 means that face j is incident upon edge i
  % so all FE(:,j) == 1 are neighbors
  FE = sparse( ...
    IC(:), ...
    repmat(1:size(tris,1),3,1)', ...
    reshape(repmat(1:3,size(tris,1),1),3*size(tris,1),1), ...
    size(sortedEdges,1), ...
    size(tris,1));

  % precompute all unit normals
  N = normalizerow(normals(vers,tris));

  adjAngles = sparse(size(tris,1),size(tris,1));
  C = sparse(size(tris,1),size(tris,1));
  
  % Peal off one set of pairs at a time
  % There is a chance this is significantly faster by first transposing FE
  while nnz(FE) > 0
    % Get index of first face per row
    [OM,J] = max(FE,[],2);
    I = 1:size(sortedEdges,1); % Edge index
    I = I(OM~=0);
    J = J(OM~=0); % First face index
    M = OM(OM~=0);
    % Lookup: L(i) = j  --> edge i reveals face j
    L = sparse(I,1,J,size(sortedEdges,1),1);
    % remove these from FE
    old = nnz(FE);
    FE = FE - sparse(I,J,M,size(sortedEdges,1),size(tris,1));
    new = nnz(FE);
    assert(new<=old);
    [I2,J2,M2] = find(FE);
    
    % Pair up with first
    J1 = L(I2);
    
    % get unit normals
    N1 = N(J1,:);
    N2 = N(J2,:);
    
    % get unit edge vector
    C1 = mod(OM(I2)+2,3)+1;
    C2 = mod(M2+2,3)+1;
    E1 = sub2ind(size(tris),J2,mod(M2+1,3)+1);
    E2 = sub2ind(size(tris),J2,mod(M2-1+1,3)+1);
    EV = normalizerow(vers(tris(E2),:)-vers(tris(E1),:));
    
    % Don't need unit normals
    % http://en.wikipedia.org/wiki/Dihedral_angle#Alternative_definitions
    D12 = pi-atan2(dot(cross(N1,N2,2),EV,2),dot(N1,N2,2));
    
    %D12 = pi-acos(sum(N1.*N2,2));
    % append to A: plus is OK here because if two facets share more than one
    % edge they are coplanar so D12 equals 0 in both cases
    adjAngles = adjAngles + sparse([J1;J2],[J2;J1],[D12;D12],size(tris,1),size(tris,1));
    C = C + sparse([J2 J1],[J1 J2],[C2 C1],size(tris,1),size(tris,1));
  end

end
