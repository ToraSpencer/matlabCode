clc;
clear all ;


%%
pathName = './test.obj';
origin = [0 0 0];
dir = [1 0 0];

length = 10;
SR = 0.5;
versCount = round(length/SR);
line = repmat(origin, versCount, 1);
addMat = SR * (0 : versCount-1)';
addMat = addMat * dir;
line = line + addMat;
OBJwriteVertices(pathName, line);

