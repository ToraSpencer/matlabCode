function printDir( pathName, dir, origin )
% 输出方向指示线点集，写入到文件。
length = 10;
SR = 0.5;
versCount = round(length/SR);
line = repmat(origin, versCount, 1);
addMat = (0 : versCount-1)';
addMat = addMat * dir;
line = line+addMat;
OBJwriteVertices(pathName, line);
end

