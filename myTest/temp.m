clc;
close all;
clear all;
%%
load('cutPatientVers.mat');
load('rootCutVers.mat');

[tris1] = MyCrustOpen(cutPatientVers);
[tris2] = MyCrustOpen(rootCutVers);

writeOBJ('�и�����.Obj', cutPatientVers, tris1);          
writeOBJ('�и�����.Obj', rootCutVers, tris2);    

%%
% ����
hole = select_holes_and_boundary(cutPatientVers, tris1);
newTris = fill_mesh_holes(cutPatientVers, tris1, hole,'closed',99999999);
newTris = double(newTris);
newTris = reduceWrongTris(newTris);  

hole = select_holes_and_boundary(cutPatientVers, newTris);
newTris = fill_mesh_holes(cutPatientVers, newTris, hole,'closed',99999999);
newTris = double(newTris);
newTris = reduceWrongTris(newTris);  


writeOBJ('�������и�����.obj', cutPatientVers, newTris);

