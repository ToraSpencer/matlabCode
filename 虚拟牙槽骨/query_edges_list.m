function [edgeList] = query_edges_list(tris, mode)
 
% ��������������ʾһ���ߡ�2�У�����Ϊ������
 
L = cat(2,tris,tris(:,1));  
R = repelem(L,1,[1 2 2 1]); % ����Ƭtris�е�[x,y,z]��дΪ[x,y,y,z,z,x];

edgeList = cell2mat(cellfun(@(x) reshape(x,[2,3])',num2cell(R,2),'un',0));

if nargin  > 1 && strcmpi(mode,'sorted')
    
    edgeList = sort(edgeList,2);   
   
elseif nargin  < 2 || strcmpi(mode,'raw')
 
    
else
    
 
end
 

end  