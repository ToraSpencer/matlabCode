clc;
close all;
clear all;
%%
load('cutPatientVers.mat');
load('rootCutVers.mat');

[tris1] = MyCrustOpen(cutPatientVers);
[tris2] = MyCrustOpen(rootCutVers);

writeOBJ('ÇÐ¸îÑÀ¹Ú.Obj', cutPatientVers, tris1);          
writeOBJ('ÇÐ¸îÑÀ¸ù.Obj', rootCutVers, tris2);    

%%
% ²¹¶´
hole = select_holes_and_boundary(cutPatientVers, tris1);
newTris = fill_mesh_holes(cutPatientVers, tris1, hole,'closed',99999999);
newTris = double(newTris);
newTris = reduceWrongTris(newTris);  

hole = select_holes_and_boundary(cutPatientVers, newTris);
newTris = fill_mesh_holes(cutPatientVers, newTris, hole,'closed',99999999);
newTris = double(newTris);
newTris = reduceWrongTris(newTris);  


writeOBJ('²¹¶´ºóÇÐ¸îÑÀ¹Ú.obj', cutPatientVers, newTris);

