function [b,bc] = boundary_conditions(V,F,C,P,E,CE)
  % BOUNDARY_CONDITIONS
  % Compute boundary and boundary conditions for solving for correspondences
  % weights over a set of mesh (V,F) with control points C(p,:) and control
  % bones C(E(:,1),:) --> C(E(:,2),:)
  %
  % [b,bc] = boundary_conditions(V,F,C,P,E,CE)
  %
  % % same as [b,bc] = boundary_conditions(V,F,C,1:size(C,1),[])
  % [b,bc] = boundary_conditions(V,F,C) 
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices (not being used...)
  %  C  list of control vertex positions 
  %  P  list of indices into C for point controls, { 1:size(C,1) }
  %  E  list of bones, pairs of indices into C, connecting control vertices, 
  %    { [] }
  %  CE  list of "cage edges", pairs of indices into ***P***, connecting
  %    control ***points***. A "cage edge" just tells point boundary conditions 
  %    to vary linearly along straight lines between the end points, and to be
  %    zero for all other handles. { [] }
  % Outputs
  % boundary_conditions.m
%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% b为边界点的索引
% bc 为矩阵，每行代表某一边界点（具体是哪一边界点由b来索引确定），每列代表某一控制            
% 点，矩阵内容代表某一控制点对某一边界点的影响权值。所有控制点对一边界点的影响权
% 值之和应该为1，所以每行按列相加也应该为1,当然所有控制点对任意顶点的影响权值之和      
% 应该为1,这里的边界点属于特殊的顶点 
%
% 注意,这里我想解释下边界点的含义,边界点其实是顶点的一种,这里用户用鼠标按下的点
% 为控制点,那么离控制点最近的点变成了我们的边界点,因为边界点必须为mesh上的顶点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  b  list of boundary vertices
  %  bc list of boundary conditions, size(boundary) by # handles matrix of
  %    boundary conditions where bc(:,i) are the boundary conditions for the 
  %    ith handle ( handle order is point handles then edges handles: P,E)
  % 
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: biharmonic_bounded 
  %

  % control vertices and domain should live in same dimensions
  assert(size(C,2) == size(V,2));
  % number o dimensions
  dim = size(C,2);

  % set point handle indices and arrange to be column vector
  if(exist('P','var'))
    if(size(P,1) == 1)
      P = P';
    end
  else
    % if point handles weren't given then treat all control vertices as point
    % handles
    P = (1:size(C,1))';
  end

  % set default edge list to []
  if(~exist('E','var'))
    E = [];
  end

  % set default cage edge list to []
  if(~exist('CE','var'))
    CE = [];
  end

  % P should be either empty or a column vector
  assert(isempty(P) || (size(P,2) == 1));
  % E should be empty or be #bones by 2 list of indices
  assert( isempty(E) || (size(E,2) == 2));


  % number of point controls
  np = numel(P);
  % number of bone controls
  ne = size(E,1);
  % number of control handles
  m = np + ne;
  % number of mesh vertices
  n = size(V, 1);

  % number of control vertices
  c = size(C,1);

  % compute distance from every vertex in the mesh to every control vertex
  D = permute(sum((repmat(V,[1,1,c]) - ...
    permute(repmat(C,[1,1,n]),[3,2,1])).^2,2),[1,3,2]);
  % 这里我想稍微解释下,假设V为n*2,那么repmat(V,[1,1,c])则为n*2*c的三维数组,
  % 其物理意义为:前两维n*2包含了n个顶点,然后把他复制成了C份
  % 那么repmat(C,[1,1,n])则为c*2*n,其物理意义为:c*2包含了c个控制点,然后复制成
  % n份,然后呢permute将其第1维和第3维置换了一下,由于第3维表示复制,现在置换过后
  % 第1维表示复制, permute(repmat(C,[1,1,n]),[3,2,1])则为n*2*c,现在前两维
  % 表示一个control point,纵向复制n份,然后控制点随着第三维变化,相减后算他们的
  % 平方,相当于每个控制点与各个顶点的距离平方,注意这里sum(*,2)是按水平方向叠加,
  %这里计算过后的维数为n*1*c,最后再进行转置,变成n*c*1,这样就变成了每个顶点到控制
  %点的距离的平方
% vrepmat
% Replicate and tile array

%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% now D is a n*c*1 matrix
%             控制点
%              －－
%      顶点|  每个点到控制点的距离平方
%          |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % use distances to determine closest mesh vertex to each control vertex
  % Cv(i) is closest vertex in V to ith control vertex in C
  [minD,Cv] = min(D);
%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% minD 为每列最小的元素
% Cv中是指哪一个顶点到该控制点距离最小,存储的是该顶点的索引
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % if number of unique closest mesh vertices is less than total number, then
  % we have contradictory boundary conditions
  if(~all(size(unique(Cv)) == size(Cv)))
    warning('Multiple control vertices snapped to the same domain vertex');
  end

  % boundary conditions for all vertices NaN means no boundary conditions
  bc = repmat(NaN,[n m]);

%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% 让在控制点的权值为1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute boundary conditions for point handles
  if( ~isempty(P) )
    bc(Cv(P),:) = eye(np,m);
  % 这里先用P当索引去引用Cv,找到边界点的索引,然后用边界点的索引去引用对应的行  
  end

%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% 让在骨骼上的点权值为1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Compute boundary conditions for bone edge handles
  if(~isempty(E))

    % average edges length to give an idea of "is on edge" tolerance
    h = avgedge(V,F);
    % loop over bones
    for( ii = 1:size(E,1) )
      [t,sqr_d] = project_to_lines(V,V(Cv(E(ii,1)),:),V(Cv(E(ii,2)),:));
      on_edge = ((abs(sqr_d) < h*1e-6) & ((t > -1e-10) & (t < (1+1e-10))));
      % get rid of any NaNs on these rows
      % WARNING: any (erroneous) point handle boundary conditions will get
      % "normalized" with bone boundary conditions
      old = bc(on_edge,:);
      old(isnan(old)) = 0;
      bc(on_edge,:) = old;
      %seamanj:这里逐渐将NaN变为0,这样就不会重复变0了,没太大的实际意义
      %bc(on_edge,:) = isnan(bc(on_edge,:)).*0 + ~isnan(bc(on_edge,:)).*bc(on_edge,:);
      bc(on_edge,np+ii) = 1;
      %seamanj:行序表示在骨骼上的顶点序号,列序先Pass掉np,即上面控制点的个数,后面是骨骼的序号
    end
  end
%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% 让在边界上的点的权值线性变化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute boundary conditions due to cage edges
  if(~isempty(CE))

    % loop over cage edges
    for( ii = 1:size(CE,1) )
      [t,sqr_d] = project_to_lines(V,V(Cv(P(CE(ii,1))),:),V(Cv(P(CE(ii,2))),:));
      h = avgedge(V,F);
      on_edge = ((abs(sqr_d) < h*1e-6) & ((t > -1e-10) & (t < (1+1e-10))));
      % get rid of any NaNs on these rows
      % WARNING: clobbers any (erroneous) point handle boundary conditions on
      % points that are on bones)
      old = bc(on_edge,:);
      old(isnan(old)) = 0;
      bc(on_edge,:) = old;
      bc(on_edge,CE(ii,1)) = 1 - t(on_edge);
      bc(on_edge,CE(ii,2)) = t(on_edge);
      %注意如果该点在边界上,那么转换成边界两端点对该点的权值,Cv(P(CE(ii,2))
      %CE代表边端点在handle里面的索引,然后用该索引去索引handle点
    end

  end

  indices = 1:n;
  % boundary is only those vertices corresponding to rows with at least one non
  % NaN entry
  b = indices(any(~isnan(bc),2));
  bc = bc(b,:);
  % replace NaNs with zeros
  bc(isnan(bc)) = 0;
%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% 归一化处理  满足所有控制点对点p的权值和为1
% A = 1     2
%     3     4
% sum(m,2) =  对每行按列求和
% 3
% 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  bc(any(bc,2),:)  = bc(any(bc,2),:) ./ repmat(sum(bc(any(bc,2),:),2),1,size(bc,2));
end