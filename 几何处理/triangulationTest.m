%% �����ʷ֣�


%% ʹ��triangle����������ʷ�
clc;
clear all;
load('./data/innerCircle.mat');
load('./data/outerCircle.mat');
load('./data/edgeInPlane.mat');
versInPlane = [innerCircle; outerCircle];
[~, trisInPlane] = triangle(versInPlane, edgeInPlane, mean(innerCircle), 'NoBoundarySteiners');
drawMesh(versInPlane, trisInPlane);
 
