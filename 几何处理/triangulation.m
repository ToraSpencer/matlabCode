%% �����ʷ֣�


%% ʹ��triangle����������ʷ�

clc;
clear all;
load('innerCircle.mat');
load('outerCircle.mat');
load('edgeInPlane.mat');
versInPlane = [innerCircle; outerCircle];
[~, trisInPlane] = triangle(versInPlane, edgeInPlane, mean(innerCircle), 'NoBoundarySteiners');
drawMesh(versInPlane, trisInPlane);
