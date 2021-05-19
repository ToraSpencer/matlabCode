function [finalVers] = biharm_solve_with_factor_modified( ...
  bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, finalVers, Omega, N0, N1)
  all = [N0 Omega];
  b0 = finalVers(N0,:);
  f1 = finalVers(N1,:);

  for coord_index = 1:3
      rhs_Dx = -bi_S(all,N0)*b0(:,coord_index) - bi_S(all,N1)*f1(:,coord_index);
      rhs_Dy = zeros(size(Omega,2),1);
      rhs = [ rhs_Dx; rhs_Dy];   
      sol = bi_Q*(bi_U\(bi_L\(bi_P*(bi_R\rhs))));
      ny = size(all,2);
      finalVers(Omega,coord_index) = sol(ny+1:end);   
  end
  
end
