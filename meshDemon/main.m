clc
close all
clear all
functionname='testforfusion.m'; 
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir 'mesh process'])


%% Read_Obj()������OBJ�ļ��ж�ȡ��������
tooth =  Read_Obj('tooth.obj');      


%% ����ά����
figure(1);
trimesh(tooth.face,tooth.vertex(:,1),tooth.vertex(:,2),tooth.vertex(:,3));
 
%% ReadObj()������OBJ�ļ��ж�ȡ��������
toothVers = ReadObj('tooth.obj');
figure(2);

%% MyCrustOpen()����������ƣ��������Ƭ����������ò�����ɵ�����Ƭ�������
trimesh(MyCrustOpen(toothVers), toothVers(:, 1), toothVers(:, 2), toothVers(:, 3));
writeOBJ('������Ƭ���ɵ�tooth����.obj', toothVers, MyCrustOpen(toothVers));

%% ������
figure(3)
plot3(toothVers(:, 1), toothVers(:, 2), toothVers(:, 3), '.');

%% writeOBJ()���������������������OBJ�ļ�
writeOBJ('newTooth.obj', tooth.vertex, tooth.face);
OBJwriteVertices('toothVers.obj', toothVers);

%% getDirLine()������ȡ����ָʾ�ߵ㼯
xline = getDirLine([1, 0, 0], [0, 0, 0]);
figure();
plot3(xline(:, 1), xline(:, 2), xline(:, 3), '*');
