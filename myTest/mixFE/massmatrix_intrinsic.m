function [ bi_M] = massmatrix_intrinsic(lalala,newTris,nvert,masstype)

    l1 = lalala(:,1); l2 = lalala(:,2); l3 = lalala(:,3);
    s = (l1 + l2 + l3)*0.5;

    dblA = 2*sqrt( s.*(s-l1).*(s-l2).*(s-l3));

    i1 = newTris(1,:); i2 = newTris(2,:); i3 = newTris(3,:); 
    
    if strcmp(masstype,'full')
        i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
        j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
        offd_v = dblA/24.;
        diag_v = dblA/12.;
        v = [offd_v,offd_v, offd_v,offd_v, offd_v,offd_v, diag_v,diag_v,diag_v];  
    elseif strcmp(masstype,'barycentric')
        i = [i1 i2 i3];
        j = [i1 i2 i3];
        diag_v = dblA/6.;
        v = [diag_v,diag_v,diag_v];
    elseif strcmp(masstype,'voronoi')
      cosines = [ ...
        (lalala(:,3).^2+lalala(:,2).^2-lalala(:,1).^2)./(2*lalala(:,2).*lalala(:,3)), ...
        (lalala(:,1).^2+lalala(:,3).^2-lalala(:,2).^2)./(2*lalala(:,1).*lalala(:,3)), ...
        (lalala(:,1).^2+lalala(:,2).^2-lalala(:,3).^2)./(2*lalala(:,1).*lalala(:,2))];
      barycentric = cosines.*lalala;
      normalized_barycentric = barycentric./[sum(barycentric')' sum(barycentric')' sum(barycentric')'];
      areas = 0.25*sqrt( ...
        (lalala(:,1) + lalala(:,2) - lalala(:,3)).* ...
        (lalala(:,1) - lalala(:,2) + lalala(:,3)).* ...
        (-lalala(:,1) + lalala(:,2) + lalala(:,3)).* ...
        (lalala(:,1) + lalala(:,2) + lalala(:,3)));
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

    else 
        error('bad mass matrix type')
    end
    bi_M = sparse(i,j,v,nvert, nvert);  
end
