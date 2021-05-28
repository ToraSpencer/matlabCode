 
function edg_b = sotr_edge(Edge_b,star)
%          center = mean(edgepoints);
%          start_pointA = v(Edge_b(1,1),:);
%          start_pointB = v(Edge_b(1,2),:);
%          C = center;
         edg_b = Edge_b(star,:);
         Edge_b(star,:) = [];
%          s = dot(n,cross(start_pointA -C , start_pointB - C));
         while ~isempty(Edge_b)
             star = edg_b(end,2);
%              A = v(star,:);
             [mm,nn] = ind2sub(size(Edge_b),find(Edge_b == star));
             edg_b = [edg_b;[star,Edge_b(mm,3-nn)]];
             Edge_b(mm,:)=[];
         end
%          edg_b = [edg_b;[edg_b(end,2),edg_b(1,1)]];
end
       