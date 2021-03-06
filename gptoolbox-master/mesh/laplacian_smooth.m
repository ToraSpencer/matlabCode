function [vers,Uall] = laplacian_smooth(vers, tris, L_method, fixedVerIdx, lambda,method, S, max_iter)

  % laplace平滑

  % [vers] = laplacian_smooth(vers,F)
  % [vers] = laplacian_smooth(vers,F,L_method,b,lambda,method,S)
  % [vers,Uall] = laplacian_smooth(vers,F,L_method,b,lambda,method,S,max_iter)

  % Inputs:
  %   vers  #V x 3 matrix of vertex coordinates
  %   tris  #F x 3  matrix of indices of triangle corners
  
  %   L_method  method for laplacian
  %      'uniform'
  %      'cotan'
  
  %   fixedVerIdx  list of indices of fixed vertices
  %   设定保持坐标不变的顶点的索引。可以设定边界不变。
  
  %   lambda  diffusion speed parameter {0.1}
  
  %   method  method to use:
  %       'implicit' (default)
  %       'explicit'
  
  %   S  scalar fields to smooth (default V)
  
  % Outputs:
  %   vers  #V x 3 list of new vertex positions
  %   Uall  #V x 3 x iters 每一次迭代的结果。
  %   

 
  %% 1. 读取参数
  n = size(vers,1);
  dim = size(vers,2);


  if(~exist('L_method','var'))
    L_method = 'cotan';
  end

  if(~exist('lambda','var')) 
    lambda = 0.1;           % 进行laplacian光顺的步长
  end

  if(~exist('fixedVerIdx','var'))
    fixedVerIdx = [];
  end

  % bulid sparse identity matrix
  I = speye(n,n);

  if(~exist('method','var'))
    method = 'implicit';
  end


  h = avgedge(vers,tris);
  if(~exist('tol','var'))
    tol = 0.001;
  end

  if(~exist('max_iter','var'))
    max_iter = 1000;
  end

  if(~exist('S','var'))
    S = vers;
  end

  
 %% 2. 计算laplacian
 
 % 平面图网格应该使用uniform laplacian
  if strcmp(L_method, 'uniform')
    A = adjacency_matrix(tris);
    L = A - diag(sum(A));               % uniform laplacian
  end

  % place for factorization and symmtery flag used by min_quad_with_fixed
  P = [];
  sym = [];

  iter = 0;
  vers = S;
  vers_prev = S;
  if nargout >= 2
    Uall = [];
  end 

  % 一般的网格应该使用余切laplacian
  if strcmp(L_method, 'cotan')
    L = cotmatrix_embedded(vers,tris);
  end

  
 %% 3. 迭代
  while( iter < max_iter && (iter == 0 || max(abs(vers(:)-vers_prev(:)))>tol*h))
    vers_prev = vers;
    
    switch method
    case 'implicit'
      Q = (I-lambda*L);
      % could prefactor Q for 'uniform' case
      for d = 1:size(S,2)
        [vers(:,d), P] = min_quad_with_fixed(Q*0.5,-vers(:,d), fixedVerIdx, S(fixedVerIdx,d), [],[],P);
      end
      
    case 'explicit'
      Q = (I+lambda*L);
      vers = Q * vers;
      % 设定保持不变的顶点：
      vers(fixedVerIdx,:) = S(fixedVerIdx, :);
      
    otherwise
      error(['' method ' is not a supported smoothing method']);
    end

    if nargout >= 2
      Uall = cat(3,Uall,vers);
    end
    
    iter = iter + 1;
  end

  
end
