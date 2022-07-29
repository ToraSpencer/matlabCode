function [ isctVer ] = planeEdgeIntersect(planeVert, planeNorm, edge, meshVers)
% 求平面和边的交点；

va = meshVers(edge(1), :);
vb = meshVers(edge(2), :);

vpDis1 = (va - planeVert)* planeNorm';
vpDis2 = (vb - planeVert)* planeNorm';

dir = vb - va;
dir = dir/norm(dir);
x = -dot(va, planeNorm)/dot(dir, planeNorm);
isctVer = va + x * dir;

end

