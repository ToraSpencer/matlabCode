function printDir( pathName, dir, origin )
% �������ָʾ�ߵ㼯��д�뵽�ļ���
length = 10;
SR = 0.5;
versCount = round(length/SR);
line = repmat(origin, versCount, 1);
addMat = (0 : versCount-1)';
addMat = addMat * dir;
line = line+addMat;
OBJwriteVertices(pathName, line);
end

