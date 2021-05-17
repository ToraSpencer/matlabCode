function [Omega,N0,N1,N2,shrinking_handle ] = ...
  layers_from_handle_modified( ...
    vertex_count, ...
    newTris, ...
    exterior)
 
  Omega = 1:vertex_count;
  Omega = Omega(~ismember(Omega,exterior));
  shrinking_handle = exterior;
  growing_interior = Omega;

  interior_faces = intersect(newTris,newTris.*ismember(newTris,growing_interior),'rows');
  handle_faces = intersect(newTris,newTris.*~ismember(newTris,growing_interior),'rows');
  H0 = setdiff(newTris,union(handle_faces,interior_faces,'rows'),'rows');
  N0 = intersect(shrinking_handle,reshape(H0,1,size(H0,1)*size(H0,2)));
  
  
  growing_interior = [growing_interior N0];
  shrinking_handle = setdiff(shrinking_handle,N0);

  interior_faces = intersect(newTris,newTris.*ismember(newTris,growing_interior),'rows');
  handle_faces = intersect(newTris,newTris.*~ismember(newTris,growing_interior),'rows');
  H1 = setdiff(newTris,union(handle_faces,interior_faces,'rows'),'rows');
  N1 = intersect(shrinking_handle,reshape(H1,1,size(H1,1)*size(H1,2)));
  growing_interior = [growing_interior N1];
  shrinking_handle = setdiff(shrinking_handle,N1);

  interior_faces = intersect(newTris,newTris.*ismember(newTris,growing_interior),'rows');
  handle_faces = intersect(newTris,newTris.*~ismember(newTris,growing_interior),'rows');
  H2 = setdiff(newTris,union(handle_faces,interior_faces,'rows'),'rows');
  N2 = intersect(shrinking_handle,reshape(H2,1,size(H2,1)*size(H2,2)));
  shrinking_handle = setdiff(shrinking_handle,N2);

 

end
