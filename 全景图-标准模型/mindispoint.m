function [n,v] = mindispoint(points,vert)
for i = 1:length(points)
    dis = sqrt((vert(:,1)-points(i,2)).^2 + (vert(:,2)-points(i,1)).^2);
    mindis(i) = min(dis);
end
n = find(mindis == min(mindis));
v = points(n,:);
end