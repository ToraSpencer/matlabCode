function [ ] = printAABB( filePath, minV, maxV )
% ����AABB��Χ������Ȼ��д�뵽�ļ��У�
arrow = maxV - minV;
xlen = arrow(1);
ylen = arrow(2);
zlen = arrow(3);

vers = [];
[vers, tris] = cube(2, 2, 2);
versT = diag([xlen, ylen, zlen]) * vers';
vers =  versT';

vers = bsxfun(@plus, vers, minV);
writeOBJ(filePath, vers, tris);


end

