function [edgeList] = query_edges_list(tris, mode)
 
% 两个顶点索引表示一条边。2列，行数为边数；
 
L = cat(2,tris,tris(:,1));  
R = repelem(L,1,[1 2 2 1]); % 三角片tris中的[x,y,z]改写为[x,y,y,z,z,x];

edgeList = cell2mat(cellfun(@(x) reshape(x,[2,3])',num2cell(R,2),'un',0));

if nargin  > 1 && strcmpi(mode,'sorted')
    
    edgeList = sort(edgeList,2);   
   
elseif nargin  < 2 || strcmpi(mode,'raw')
 
    
else
    
 
end
 

end  