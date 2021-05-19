function bi_M = massmatrix(mergedToothVers,newTris, type)

    i1 = newTris(1,:); 
    i2 = newTris(2,:); 
    i3 = newTris(3,:); 
    v1 = mergedToothVers(i3,:) - mergedToothVers(i2,:);  
    v2 = mergedToothVers(i1,:) - mergedToothVers(i3,:); 
    v3 = mergedToothVers(i2,:) - mergedToothVers(i1,:);

    if size(mergedToothVers,2) == 2
      dblA = v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1);
    elseif size(mergedToothVers,2) == 3
      n  = cross(v1,v2,2);  

      dblA = (sqrt(sum((n').^2)))';
    else 
      error('unsupported vertex dimension %d', size(mergedToothVers,2))
    end
    if strcmp(type,'full')
        i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
        j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
        offd_v = dblA/24.;
        diag_v = dblA/12.;
        v = [offd_v,offd_v, offd_v,offd_v, offd_v,offd_v, diag_v,diag_v,diag_v];  
        bi_M = sparse(i,j,v,size(mergedToothVers,1), size(mergedToothVers,1));
    elseif strcmp(type,'barycentric')
        % only diagonal elements
        i = [i1 i2 i3];
        j = [i1 i2 i3];
        diag_v = dblA/6.;
        v = [diag_v,diag_v,diag_v];
        bi_M = sparse(i,j,v,size(mergedToothVers,1), size(mergedToothVers,1));
    elseif strcmp(type,'voronoi')

      FT = newTris';
      l = [ ...
        sqrt(sum((mergedToothVers(FT(:,2),:)-mergedToothVers(FT(:,3),:)).^2,2)) ...
        sqrt(sum((mergedToothVers(FT(:,3),:)-mergedToothVers(FT(:,1),:)).^2,2)) ...
        sqrt(sum((mergedToothVers(FT(:,1),:)-mergedToothVers(FT(:,2),:)).^2,2)) ...
        ];
      bi_M = massmatrix_intrinsic(l,newTris,size(mergedToothVers,1),'voronoi');
    else 
        error('bad mass matrix type')
    end
end
