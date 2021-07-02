%% 三角剖分：


%% 使用triangle库进行三角剖分
clc;
clear all;
load('innerCircle.mat');
load('outerCircle.mat');
load('edgeInPlane.mat');
versInPlane = [innerCircle; outerCircle];
[~, trisInPlane] = triangle(versInPlane, edgeInPlane, mean(innerCircle), 'NoBoundarySteiners');
drawMesh(versInPlane, trisInPlane);


%% 使用triangle库进行三角剖分
clc;
clear all;
vers = readOBJ('twoCircles.obj');
versCount = size(vers, 1);
versInPlane = vers(:, 1:2);
edgeInPlane = zeros(versCount, 2);
edgeInPlane(:,1) = 1:versCount;
edgeInPlane(1:29,2) = 2:30;
edgeInPlane(30, 2) = 1;
edgeInPlane(31:59, 2) = 32:60;
edgeInPlane(60, 2) = 31;
[~, trisInPlane] = triangle(versInPlane, edgeInPlane, mean(versInPlane), 'NoBoundarySteiners');
drawMesh(versInPlane, trisInPlane);
