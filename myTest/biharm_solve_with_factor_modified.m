function [finalVers] = biharm_solve_with_factor_modified( ...
                A, bi_S, finalVers, Omega, N0, N1)
  all = [N0 Omega];
  b0 = finalVers(N0,:);
  f1 = finalVers(N1,:);

  for coord_index = 1:3
      rhs_Dx = -bi_S(all,N0)*b0(:,coord_index) - bi_S(all,N1)*f1(:,coord_index);
      rhs_Dy = zeros(size(Omega,2),1);
      rhs = [ rhs_Dx; rhs_Dy];   
      
      x = A\rhs;
    % sol = bi_Q*( inv(bi_U) * ( inv(bi_L) * ( bi_P * ( inv(bi_D) * rhs))));
    % D * inv(bi_P) * bi_L * bi_U * inv(bi_Q) * sol == rhs;
    % A * x == rhs;
    % x是上面线性方程组的解，是一条列向量
    % 原程序使用LU分解的方法来解线性方程组，这里直接使用左除运算。
      
      ny = size(all,2);
      finalVers(Omega,coord_index) = x(ny+1:end);   
  end
  
  
  disp('finished.');
end
