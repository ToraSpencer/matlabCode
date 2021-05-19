function bi_S = cotmatrix(mergedToothVers,newTris)
  i1 = newTris(1,:); i2 = newTris(2,:); i3 = newTris(3,:); 
  v1 = mergedToothVers(i3,:) - mergedToothVers(i2,:);  v2 = mergedToothVers(i1,:) - mergedToothVers(i3,:); v3 = mergedToothVers(i2,:) - mergedToothVers(i1,:);
  n  = cross(v1,v2,2); 
  dblA = (sqrt(sum((n').^2)))';
  cot12 = -dot(v1,v2,2)./dblA/2; cot23 = -dot(v2,v3,2)./dblA/2; cot31 = -dot(v3,v1,2)./dblA/2;
  diag1 = -cot12-cot31; diag2 = -cot12-cot23; diag3 = -cot31-cot23;
  i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
  j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
  v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];
  bi_S = sparse(i,j,v,size(mergedToothVers,1),size(mergedToothVers,1));
end
