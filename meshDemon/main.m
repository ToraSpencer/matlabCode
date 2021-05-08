clc
close all
clear all
functionname='testforfusion.m'; 
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir 'mesh process'])


%% Read_Obj()——从OBJ文件中读取网格数据
tooth =  Read_Obj('tooth.obj');      


%% 画三维网格
figure(1);
trimesh(tooth.face,tooth.vertex(:,1),tooth.vertex(:,2),tooth.vertex(:,3));
 
%% ReadObj()——从OBJ文件中读取点云数据
toothVers = ReadObj('tooth.obj');
figure(2);

%% MyCrustOpen()——输出点云，输出三角片——！！！貌似生成的三角片朝向混乱
trimesh(MyCrustOpen(toothVers), toothVers(:, 1), toothVers(:, 2), toothVers(:, 3));
writeOBJ('补三角片生成的tooth网格.obj', toothVers, MyCrustOpen(toothVers));

%% 画点云
figure(3)
plot3(toothVers(:, 1), toothVers(:, 2), toothVers(:, 3), '.');

%% writeOBJ()——网格或点云数据输出到OBJ文件
writeOBJ('newTooth.obj', tooth.vertex, tooth.face);
OBJwriteVertices('toothVers.obj', toothVers);

%% getDirLine()——获取方向指示线点集
xline = getDirLine([1, 0, 0], [0, 0, 0]);
figure();
plot3(xline(:, 1), xline(:, 2), xline(:, 3), '*');
