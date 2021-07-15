function [ B ] = normalizerow( A ) %#codegen
  % NORMALIZEROW Normalize each row so that each row's l2 norm as a vector is 1
  %
  % [ B ] = normalizerow( A )
  %
  % Input:
  %  A  #A by D list of row vectors of dimension D
  % Output:
  %  B  #B by D list of normalized row vectors 
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), Daniele Panozzo
  %

  if issparse(A)
    % speed up (20x) for large sparse matrices
    B = bsxfun(@times,A,1./sqrt(sum(A.^2,2)));
  else
    %% my debug
    %      B = A./normrow(A);
    temp = normrow(A);                      % ������
    temp = ones(size(temp))./temp;
    temp  = repmat(temp, 1, size(A, 2));
    B = A.*temp;
    
  end
end

