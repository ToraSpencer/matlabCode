function [ normDirs ] = normals(vers, tris, varargin)
    % 计算多边形网格的面片的法向（未归一化）；
    
  % NORMALS Compute *unnormalized* normals per face
  % N = normals(V,F)
  % N = normals(V,F,'ParameterName',ParameterValue, ...)
 
  % Inputs:
  %  V  #V x 3 matrix of vertex coordinates
  %  F  #F x 3  matrix of indices of triangle corners
  %  Optional:
  %    'Stable' followed by whether to compute normals in a way stable with
  %      respect to vertex order: constant factor more expensive {false}
  %    'UseSVD' followed by whether to use SVD, slow {false}
  % Output:
  %  N  #F x 3 list of face normals
 
  % Example:
  %   timeit(@() normalizerow(normals(UV,UF)))
  %   % faster unit normals
  %   cross2 = @(a,b,c) ...
  %     [a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
  %      a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
  %      a(:,1).*b(:,2)-a(:,2).*b(:,1)];
  %   nrmlev = @(p1,p2,p3) cross2(p2-p1,p3-p1);
  %   unrml = @(V,F) normalizerow(nrmlev(V(F(:,1),:),V(F(:,2),:),V(F(:,3),:)));
  %   timeit(@() unrml(UV,UF))


  function D = sum3(A,B,C)
    % SUM3 Entrywise sum of three matrices in a stable way: sorting entries by
    % value and summing 
    shape = size(A);
    ABC = [A(:) B(:) C(:)];
    [~,I] = sort(abs(ABC),2,'descend');
    sABC = reshape( ...
      ABC(sub2ind(size(ABC),repmat(1:size(ABC,1),size(ABC,2),1)',I)),size(ABC));
    D = reshape(sum(sABC,2),shape);
  end

  % 默认optional参数：
  stable = false;
  use_svd = false;
  
  % 读取optional参数；
  % default values, Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Stable','UseSVD'},{'stable','use_svd'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables, param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller', params_to_variables(param_name), varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  
  if size(tris, 2) == 2
    normDirs = vers(tris(:,2),:)-vers(tris(:,1),:);
    normDirs = normDirs*[0 -1;1 0];
    return;
  end
  
  va = vers(tris(:,1),:);
  vb = vers(tris(:,2),:);
  vc = vers(tris(:,3),:);
  
  if use_svd                       % 'UseSVD', 1——使用SVD求面片法向；
    normDirs = zeros(size(tris, 1),3);
    BC = barycenter(vers, tris);
    for k = 1:size(tris,1)
      Uf = bsxfun(@minus, vers(tris(k,:),:), BC(k,:));
      [~,~,sV] = svd(Uf);
      normDirs(k,:) = sV(:,3);
    end
    NN = normals(vers,tris,'UseSVD',false);
    normDirs(sum(normDirs.*NN,2)<0,:) = normDirs(sum(normDirs.*NN,2)<0,:)*-1;
  else
    % 无optional参数
    N1 = cross(vb - va, vc - va, 2);        % 2 is necessary because this will produce the wrong result if there are exactly 3 faces
    
    if stable         % 'Stable', 1
      N2 = cross(vc - vb, va - vb, 2);
      N3 = cross(va - vc, vb - vc, 2);
      normDirs = sum3(N1,N2,N3)/3;
    else
      normDirs = N1;
    end
  end


end

