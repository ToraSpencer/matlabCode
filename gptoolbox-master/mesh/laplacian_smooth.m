function [U,Uall] = laplacian_smooth(vers, tris, L_method, b, lambda,method, S, max_iter)

  % laplace平滑

  % [U] = laplacian_smooth(V,F)
  % [U] = laplacian_smooth(V,F,L_method,b,lambda,method,S)
  % [U,Uall] = laplacian_smooth(V,F,L_method,b,lambda,method,S,max_iter)

  % Inputs:
  %   V  #V x 3 matrix of vertex coordinates
  %   F  #F x 3  matrix of indices of triangle corners
  
  %   L_method  method for laplacian
  %      'uniform'
  %      'cotan'
  
  %   b  list of indices of fixed vertices ？？？ 用例里面传入的是空向量。
  
  %   lambda  diffusion speed parameter {0.1}
  
  %   method  method to use:
  %       'implicit' (default)
  %       'explicit'
  
  %   S  scalar fields to smooth (default V)
  
  % Outputs:
  %   U  #V x 3 list of new vertex positions
  %   Uall  #V x 3 x iters 每一次迭代的结果。
  %   

 
  %% 1. 读取参数
  n = size(vers,1);
  dim = size(vers,2);


  if(~exist('L_method','var'))
    L_method = 'cotan';
  end

  if(~exist('lambda','var')) 
    lambda = 0.1; 
  end

  if(~exist('b','var'))
    b = [];
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
  U = S;
  U_prev = S;
  if nargout >= 2
    Uall = [];
  end 

  % 余切laplacian
  if strcmp(L_method,'cotan')
    L = cotmatrix_embedded(vers,tris);
  end

  
 %% 3. 迭代
  while( iter < max_iter && (iter == 0 || max(abs(U(:)-U_prev(:)))>tol*h))
    U_prev = U;
    switch method
    case 'implicit'
      Q = (I-lambda*L);
      % could prefactor Q for 'uniform' case
      for d = 1:size(S,2)
        [U(:,d),P] = min_quad_with_fixed(Q*0.5,-U(:,d),b,S(b,d),[],[],P);
      end
    case 'explicit'
      Q = (I+lambda*L);
      U = Q * U;
      % enforce boundary
      U(b,:) = S(b,:);
    otherwise
      error(['' method ' is not a supported smoothing method']);
    end

    if nargout >= 2
      Uall = cat(3,Uall,U);
    end
    iter = iter + 1;
  end

  
end
