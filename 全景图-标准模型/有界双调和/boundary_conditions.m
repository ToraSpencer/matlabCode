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
% bΪ�߽�������
% bc Ϊ����ÿ�д���ĳһ�߽�㣨��������һ�߽����b������ȷ������ÿ�д���ĳһ����            
% �㣬�������ݴ���ĳһ���Ƶ��ĳһ�߽���Ӱ��Ȩֵ�����п��Ƶ��һ�߽���Ӱ��Ȩ
% ֵ֮��Ӧ��Ϊ1������ÿ�а������ҲӦ��Ϊ1,��Ȼ���п��Ƶ�����ⶥ���Ӱ��Ȩֵ֮��      
% Ӧ��Ϊ1,����ı߽����������Ķ��� 
%
% ע��,������������±߽��ĺ���,�߽����ʵ�Ƕ����һ��,�����û�����갴�µĵ�
% Ϊ���Ƶ�,��ô����Ƶ�����ĵ��������ǵı߽��,��Ϊ�߽�����Ϊmesh�ϵĶ���
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
  % ����������΢������,����VΪn*2,��ôrepmat(V,[1,1,c])��Ϊn*2*c����ά����,
  % ����������Ϊ:ǰ��άn*2������n������,Ȼ��������Ƴ���C��
  % ��ôrepmat(C,[1,1,n])��Ϊc*2*n,����������Ϊ:c*2������c�����Ƶ�,Ȼ���Ƴ�
  % n��,Ȼ����permute�����1ά�͵�3ά�û���һ��,���ڵ�3ά��ʾ����,�����û�����
  % ��1ά��ʾ����, permute(repmat(C,[1,1,n]),[3,2,1])��Ϊn*2*c,����ǰ��ά
  % ��ʾһ��control point,������n��,Ȼ����Ƶ����ŵ���ά�仯,����������ǵ�
  % ƽ��,�൱��ÿ�����Ƶ����������ľ���ƽ��,ע������sum(*,2)�ǰ�ˮƽ�������,
  %�����������ά��Ϊn*1*c,����ٽ���ת��,���n*c*1,�����ͱ����ÿ�����㵽����
  %��ľ����ƽ��
% vrepmat
% Replicate and tile array

%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% now D is a n*c*1 matrix
%             ���Ƶ�
%              ����
%      ����|  ÿ���㵽���Ƶ�ľ���ƽ��
%          |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % use distances to determine closest mesh vertex to each control vertex
  % Cv(i) is closest vertex in V to ith control vertex in C
  [minD,Cv] = min(D);
%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% minD Ϊÿ����С��Ԫ��
% Cv����ָ��һ�����㵽�ÿ��Ƶ������С,�洢���Ǹö��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % if number of unique closest mesh vertices is less than total number, then
  % we have contradictory boundary conditions
  if(~all(size(unique(Cv)) == size(Cv)))
    warning('Multiple control vertices snapped to the same domain vertex');
  end

  % boundary conditions for all vertices NaN means no boundary conditions
  bc = repmat(NaN,[n m]);

%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% ���ڿ��Ƶ��ȨֵΪ1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute boundary conditions for point handles
  if( ~isempty(P) )
    bc(Cv(P),:) = eye(np,m);
  % ��������P������ȥ����Cv,�ҵ��߽�������,Ȼ���ñ߽�������ȥ���ö�Ӧ����  
  end

%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% ���ڹ����ϵĵ�ȨֵΪ1
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
      %seamanj:�����𽥽�NaN��Ϊ0,�����Ͳ����ظ���0��,û̫���ʵ������
      %bc(on_edge,:) = isnan(bc(on_edge,:)).*0 + ~isnan(bc(on_edge,:)).*bc(on_edge,:);
      bc(on_edge,np+ii) = 1;
      %seamanj:�����ʾ�ڹ����ϵĶ������,������Pass��np,��������Ƶ�ĸ���,�����ǹ��������
    end
  end
%%%%%%%% Added by seamanj %%%%%%%%%%%%%%
% ���ڱ߽��ϵĵ��Ȩֵ���Ա仯
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
      %ע������õ��ڱ߽���,��ôת���ɱ߽����˵�Ըõ��Ȩֵ,Cv(P(CE(ii,2))
      %CE����߶˵���handle���������,Ȼ���ø�����ȥ����handle��
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
% ��һ������  �������п��Ƶ�Ե�p��Ȩֵ��Ϊ1
% A = 1     2
%     3     4
% sum(m,2) =  ��ÿ�а������
% 3
% 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  bc(any(bc,2),:)  = bc(any(bc,2),:) ./ repmat(sum(bc(any(bc,2),:),2),1,size(bc,2));
end