function bdryLength = boundary_length(vers, tris)
% �����Ե���߳���
 
bdryEdges = outline(tris);

 
bdryEdgesLens = vers(bdryEdges(:,2),:) - vers(bdryEdges(:,1),:);

%  temp =  normrow(bdryEdgesLens);
temp = bdryEdgesLens.^2;
temp = sum(temp, 2);
temp = sqrt(temp);
bdryLength = sum(temp);

end

