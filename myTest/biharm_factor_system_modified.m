function [bi_L,bi_U,bi_P,bi_Q,bi_R, bi_S, bi_M] = biharm_factor_system_modified( ...
  mergedToothVers, ...
  newTris, ...
  omega, ...
  N0, ...
  N1)

  interest = [omega N0 N1 ];
  
  newTris = newTris( ...
    ismember(newTris(:,1),interest) & ...
    ismember(newTris(:,2),interest) & ...
    ismember(newTris(:,3),interest),:)';  

  all = [N0 omega];
  
  %bi_S = cotmatrix(mergedToothVers, newTris);  
  i1 = newTris(1,:); i2 = newTris(2,:); i3 = newTris(3,:); 
  vers1 = mergedToothVers(i3,:) - mergedToothVers(i2,:);  
  vers2 = mergedToothVers(i1,:) - mergedToothVers(i3,:); 
  vers3 = mergedToothVers(i2,:) - mergedToothVers(i1,:);
  vers4  = cross(vers1,vers2,2); 
  dblA = (sqrt(sum((vers4').^2)))';
  
  cot12 = -dot(vers1,vers2,2)./dblA/2; cot23 = -dot(vers2,vers3,2)./dblA/2; cot31 = -dot(vers3,vers1,2)./dblA/2;
  diag1 = -cot12-cot31; diag2 = -cot12-cot23; diag3 = -cot31-cot23;
  
  i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
  j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
  v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];
  [m,n] = size(v);
  v2 = reshape(v, m*n, 1);
  bi_S = sparse(i,j,v2,size(mergedToothVers,1),size(mergedToothVers,1));
 
 
  %bi_M = massmatrix(mergedToothVers, newTris, 'voronoi');
  elems = nonzeros(bi_S);
  theElem = bi_S(i(1), j(1));
  
  ml = [ ...
    sqrt(sum((mergedToothVers(newTris(2,:),:)-mergedToothVers(newTris(3,:),:)).^2,2)) ...
    sqrt(sum((mergedToothVers(newTris(3,:),:)-mergedToothVers(newTris(1,:),:)).^2,2)) ...
    sqrt(sum((mergedToothVers(newTris(1,:),:)-mergedToothVers(newTris(2,:),:)).^2,2)) ...
    ];
 

 %bi_M = massmatrix_intrinsic(lalala, newTris, size(mergedToothVers,1),'voronoi');
    
    i1 = newTris(1,:); i2 = newTris(2,:); i3 = newTris(3,:); 
    
      cosines = [ ...
        (ml(:,3).^2+ml(:,2).^2-ml(:,1).^2)./(2*ml(:,2).*ml(:,3)), ...
        (ml(:,1).^2+ml(:,3).^2-ml(:,2).^2)./(2*ml(:,1).*ml(:,3)), ...
        (ml(:,1).^2+ml(:,2).^2-ml(:,3).^2)./(2*ml(:,1).*ml(:,2))];
      barycentric = cosines.*ml;
      normalized_barycentric = barycentric./[sum(barycentric')' sum(barycentric')' sum(barycentric')'];
      
      
      areas = 0.25*sqrt( ...
        (ml(:,1) + ml(:,2) - ml(:,3)).* ...
        (ml(:,1) - ml(:,2) + ml(:,3)).* ...
        (-ml(:,1) + ml(:,2) + ml(:,3)).* ...
        (ml(:,1) + ml(:,2) + ml(:,3)));
    
      partial_triangle_areas = normalized_barycentric.*[areas areas areas];
      
      quads = [ (partial_triangle_areas(:,2)+ partial_triangle_areas(:,3))*0.5 ...
        (partial_triangle_areas(:,1)+ partial_triangle_areas(:,3))*0.5 ...
        (partial_triangle_areas(:,1)+ partial_triangle_areas(:,2))*0.5];
    
      quads(cosines(:,1)<0,:) = [areas(cosines(:,1)<0,:)*0.5, ...
        areas(cosines(:,1)<0,:)*0.25, areas(cosines(:,1)<0,:)*0.25];
    
      quads(cosines(:,2)<0,:) = [areas(cosines(:,2)<0,:)*0.25, ...
        areas(cosines(:,2)<0,:)*0.5, areas(cosines(:,2)<0,:)*0.25];
      quads(cosines(:,3)<0,:) = [areas(cosines(:,3)<0,:)*0.25, ...
        areas(cosines(:,3)<0,:)*0.25, areas(cosines(:,3)<0,:)*0.5];
    
      i = [i1 i2 i3];
      j = [i1 i2 i3];
      v = reshape(quads,size(quads,1)*3,1);
    bi_M = sparse(i,j,v,size(mergedToothVers,1), size(mergedToothVers,1));
    
    

    
  n_Omega = size(omega,2);
  Z_Omega_Omega = sparse(n_Omega, n_Omega);
  
  temp1 = -bi_M(all,all) ;
  temp2 = bi_S(all,omega);
  temp3 = bi_S(omega,all);
  temp4 = Z_Omega_Omega;
  
  % for debug
  elems1 = nonzeros(temp1);
  elems2 = nonzeros(temp2);
  elems3 = nonzeros(temp3);
    A = [ -bi_M(all,all)      bi_S(  all,omega);   ...        % 546 X 546
         bi_S(omega,all)    Z_Omega_Omega ];     
     
  elemsa = nonzeros(A);   
     
  [bi_L,bi_U,bi_P,bi_Q,bi_R] = lu(A);
  
  disp('finished.');
end
