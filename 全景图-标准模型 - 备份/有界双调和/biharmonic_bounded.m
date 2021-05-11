function W = biharmonic_bounded(V,F,b,bc,type,pou,low,up)
  % BIHARMONIC_BOUNDED Compute biharmonic bounded coordinates, using quadratic
  % optimizer
  %
  % W = biharmonic_bounded(V,F,b,bc,type,pou)
  % W = biharmonic_bounded(V,F,b,bc,type,pou,low,up)
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices, for 3D F is #F by 4, for 2D F is #F by 3
  %seamanj added: for 3D it's tetrahedron, for 2D it's triangle%
  %  b  list of boundary vertices
  %  bc list of boundary conditions, size(boundary) by # handles matrix of
  %    boundary conditions where bc(:,i) are the boundary conditions for the 
  %    ith handle
  % 这里只处理了控制点的情况,骨骼和cage没有在此处理
  %  Optional:
  %    type  type of optimizer to use {best available}:
  %      'quad'
  %      'least-squares'
  %      'conic'
  %    pou  true or false, enforce partition of unity explicitly {false}
  %    low  lower bound {0}
  %    up  upper bound {1}
  %  
  %
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: boundary_conditions
  %

  % set default for enforcing partition of unity constraint
  if ~exist('pou','var') || isempty(pou)
    pou = false;
  end

  % number of vertices
  n = size(V,1);
  % number of handles
  m = size(bc,2);

  % Build discrete laplacian and mass matrices used by all handles' solves
  if(size(F,2)==4)
    fprintf('Solving over volume...\n');
    L = cotmatrix3(V,F);
    M = massmatrix3(V,F,'barycentric');
  else
    L = cotmatrix(V,F);
    M = massmatrix(V,F,'voronoi');
  end

  % default bounds
  if ~exist('low','var') || isempty(low)
    low = 0;
  end
  if ~exist('up','var') || isempty(up)
    up = 1;
  end

  % set default optimization method
  if ~exist('type','var') || isempty(type)
    if(exist('mosekopt'))
      % if mosek is available then conic is fastest
      type = 'conic';
    else
      % if we only have matlab then quadratic is default
      type = 'quad';
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SET UP SOLVER
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % check for mosek and set its parameters
  param = [];
  % Tolerance parameter
  % >1e0 NONSOLUTION
  % 1e-1 artifacts in deformation
  % 1e-3 artifacts in isolines
  % 1e-4 seems safe for good looking deformations
  % 1e-8 MOSEK DEFAULT SOLUTION
  % 1e-14 smallest allowed value
  if(exist('mosekopt','file'))
    if(strcmp(type,'quad') || strcmp(type,'least-squares'))
      param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-10;
    elseif(strcmp(type,'conic'))
      param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-10;
    end
    mosek_exists = true;
  else 
    mosek_exists = false;
    if(verLessThan('matlab','7.12'))
      % old matlab does not solve quadprog with sparse matrices: SLOW
      % solution: dowloand MOSEK or upgrade to 2011a or greater
      warning([ ...
        'You are using an old version of MATLAB that does not support ' ...
        'solving large, sparse quadratic programming problems. The ' ...
        'optimization will be VERY SLOW and the results will be ' ...
        'INACCURATE. Please install Mosek or upgrade to MATLAB version >= ' ...
        '2011a.']);
    else
      % Tell matlab to use interior point solver, and set tolerance
      % 1e-8 MATLAB DEFAULT SOLUTION (very low accuracy)
      % 1e-10 (low accuracy)
      % 1e-12 (medium-low accuracy)
      % 1e-14 (medium accuracy)
      % 1e-16 (high accuracy)
      param = optimset( ...
        'TolFun',1e-16, ...
        'Algorithm','interior-point-convex', ...
        ... % 'Algorithm','active-set', ...
        'MaxIter', 1000, ...
        'Display','off');
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SET UP PROBLEM AND SOLVE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(pou)
    % Enforce partition of unity as explicity constraints: solve for weights
    % of all handles simultaneously
    % seamanj: Enforce partition of unity 这里的归一性,代表所有控制点对某顶点的权值和应该为1
    if(strcmp(type,'quad'))
      % biharmonic system matrix
      Qi = L*(M\L);
      % x = A\B 
      % divides the Galois array A into B to produce a particular solution 
      % of the linear equation A*x = B. In the special case when A is a 
      % nonsingular square matrix, x is the unique solution, inv(A)*B, 
      % to the equation.
      Q = sparse(m*n,m*n);
      % Q is sparse matrix with Qi along diagonal
      for ii = 1:m
        d = (ii - 1)*n + 1;
        Q(d:(d + n-1), d:(d + n-1)) = Qi;
      end
      % linear constraints: partition of unity constraints and boundary
      % conditions
      PA = repmat(speye(n,n),1,m);
      Pb = ones(n,1);
      %      1             2                          m
      %----------------------------------------------------
      % 1             1                           1
      %   1             1                           1
      %     ...           ...           ...           ...
      %         1             1                           1
      %            1             1                           1
      %然后未知量为:
      %  w_11
      %  w_21
      %  ...
      %  w_n1
      %  w_12
      %  w_22
      %  ...
      %  w_n2
      %  ...
      %  w_1m
      %  w_2m
      %  ...
      %  w_nm
      %其中w_ij表示第j个控制点对第i个顶点的权值
      %与上面矩阵相乘则得到:w_11 + w_12 + ... + w_1m,则表示所有控制点对某一顶点的权值和应该为1
      % boundary conditions
      BCAi = speye(n,n);
      BCAi = BCAi(b,:);%b为控制点的序号,这里把控制点的行选出,这里行代表控制点,列代表顶点
      BCA = sparse(m*size(BCAi,1),m*size(BCAi,2));
      % BCA is sparse matrix with BCAi along diagonal
      for ii = 1:m
        di = (ii - 1)*size(BCAi,1) + 1;
        dj = (ii - 1)*size(BCAi,2) + 1;
        BCA(di:(di + size(BCAi,1)-1), dj:(dj + size(BCAi,2)-1)) = BCAi;
      end
      % 然后BCA将BCAi从行列扩展m倍,然后对角线矩阵为BCAi,相当于BCAi沿对角线的方向
      % 扩展了m倍
      BCb = bc(:);%每列往上一列下面叠加,形成一个列向量
      %BCA*x=BCb,x解出来的值为(n*m,1)的列向量,其中数据以列叠加在一起,每列表示
      %某一控制点对所以顶点的权值的分布
      %seamanj:
%     A = [1 2; 3 4];
%     A(:)
%     ans =
%          1
%          3
%          2
%          4
      % set bounds
      ux = up.*ones(m*n,1);
      lx = low.*ones(m*n,1);
      if(mosek_exists)
        fprintf('Quadratic optimization using mosek...\n');
      else
        fprintf('Quadratic optimization using matlab...\n');
      end
      fprintf( [ ...
        '  minimize:     x''LM\\Lx\n' ...
        'subject to: %g <= x <= %g, ∑_i xi = 1\n'], ...
        low,up);
      tic;
      W = quadprog(Q,zeros(n*m,1),[],[],[PA;BCA],[Pb;BCb],lx,ux,[],param);

      %x解出来的值为(n*m,1)的列向量,其中数据以列叠加在一起,每列表示
      %某一控制点对所以顶点的权值的分布
      toc
      W = reshape(W,n,m);
%     seamanj:      
%     A = [1:12];
%     reshape(A,3,4)
% 
%     ans =
% 
%     1     4     7    10
%     2     5     8    11
%     3     6     9    12
    else
      error( [ ...
        'Enforcing partition of unity only support in conjunction with ' ...
        'type=''quad''']);
    end
    %后面的情况差不多,我就不一一分析了,只不过用了其他的优化方式
  else 
    % Drop partition of unity constraints, solve for weights of each handle
    % independently then normalize to enforce partition of unity
    % seamanj:这里没有将归一性作为constrain,而是最后去单位化
    if(strcmp(type,'quad'))
      % build quadratic coefficient matrix (bilaplacian operator)
      Q = L*(M\L);
      % set bounds
      ux = up.*ones(n,1);
      lx = low.*ones(n,1);
    elseif(strcmp(type,'least-squares'))
      % solve same problem but as least-squares problem see mosek documention
      % for details
      I = speye(n);
      Z = sparse(n,n);
      Q = [Z,Z;Z,I];
      F = sqrt(M)\L;
      c = zeros(n,1);
      B = [F,-I];
      ux = [up.*ones(n,1) ;  Inf*ones(n,1)];
      lx = [low.*ones(n,1); -Inf*ones(n,1)];
    elseif(strcmp(type,'conic'))
      % solve same problem but as conic problem see mosek documention for
      % details
      F = sqrt(M)\L;
      prob.c = [zeros(2*n,1); 1];
      I = speye(n);
      prob.a = [F,-I,zeros(n,1)];
      prob.blc = zeros(n,1);
      prob.buc = zeros(n,1);
      prob.bux = [ up.*ones(n,1);  Inf*ones(n,1);  Inf];
      prob.blx = [ low.*ones(n,1); -Inf*ones(n,1); -Inf];
      prob.cones = cell(1,1);
      prob.cones{1}.type = 'MSK_CT_QUAD';
      t_index = 2*n +1;
      z_indices = (n+1):(2*n);
      prob.cones{1}.sub = [t_index z_indices];
    else
      error('Bad type');
    end

    % number of handles
    m = size(bc,2);
    % allocate space for weights
    W = zeros(n,m);
    tic;
    % loop over handles
    for i = 1:m
      if(strcmp(type,'quad'))
        % enforce boundary conditions via lower and upper bounds
        %lx(b) = bc(:,i);
        %ux(b) = bc(:,i);
        Aeq = speye(n,n);
        Aeq = Aeq(b,:);
        if(mosek_exists)
          fprintf('Quadratic optimization using mosek...\n');
        else  
          fprintf('Quadratic optimization using matlab...\n');
        end
        fprintf( [ ...
          '  minimize:     x''LM\\Lx\n' ...
          'subject to: %g <= x <= %g\n' ], ...
          low,up);
        % if mosek is not available, then matlab will complain that sparse
        % matrices are not yet supported...
        [x,fval,err] = quadprog(Q,zeros(n,1),[],[],Aeq,bc(:,i),lx,ux,[],param);
        if(err ~= 1)
          fprintf([...
            '----------------------------------------------------------\n' ...
            'ERROR ('  num2str(err) ',' num2str(fval) '):' ...
            ' solution may be inaccurate...\n' ...
            '----------------------------------------------------------\n' ...
            ]);
        end
      elseif(strcmp(type,'least-squares'))
        % enforce boundary conditions via lower and upper bounds
        lx(b) = bc(:,i);
        ux(b) = bc(:,i);
        fprintf('Quadratic optimization using mosek...\n');
        fprintf([ ...
          '  minimize:       z''z\n' ...
          '  subject to: sqrt(M)\\Lx - z = 0\n' ...
          '  and          %g <= x <= %g\n'] , ...
          low,up);
        x = quadprog(Q,zeros(2*n,1),[],[],B,c,lx,ux,[],param);
      elseif(strcmp(type,'conic'))
        prob.bux(b) = bc(:,i);
        prob.blx(b) = bc(:,i);
        fprintf('Conic optimization using mosek...\n');
        fprintf([ ...
          '  minimize:         t\n' ...
          '  subject to: sqrt(M)\\Lx - z = 0,\n' ...
          '             t >= sqrt(z''z),\n' ...
          '               %f <= x <= %f\n'], ...
          low,up);
        [r,res]=mosekopt('minimize echo(0)',prob,param);
        % check for mosek error
        if(r == 4006)
          warning(['MOSEK ERROR. rcode: ' ...
            num2str(res.rcode) ' ' ...
            res.rcodestr ' ' ...
            res.rmsg ...
            'The solution is probably OK, but ' ...
            'to make this error go away, increase: ' ...
            'MSK_DPAR_INTPNT_CO_TOL_REL_GAP' ...
            n]);
        elseif(r ~= 0)
          error(['FATAL MOSEK ERROR. rcode: ' ...
            num2str(res.rcode) ' ' ...
            res.rcodestr ' ' ...
            res.rmsg]);
        end
        % extract solution from result
        x = res.sol.itr.xx;
      end
      % set weights to solution in weight matrix
      W(:,i) = x(1:n);
      fprintf('Lap time: %gs\n',toc);
    end
    t = toc;
    fprintf('Total elapsed time: %gs\n',t);
    fprintf('Average time per handle: %gs\n',t/m);
  end
end