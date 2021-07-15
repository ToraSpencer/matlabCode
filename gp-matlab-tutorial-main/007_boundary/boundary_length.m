function bdryLength = boundary_length(vers, tris)
% 计算边缘曲线长度
 
bdryEdges = outline(tris);

 
bdryEdgesLens = vers(bdryEdges(:,2),:) - vers(bdryEdges(:,1),:);

%  temp =  normrow(bdryEdgesLens);
temp = bdryEdgesLens.^2;
temp = sum(temp, 2);
temp = sqrt(temp);
bdryLength = sum(temp);

end

