%% 三角剖分：


%% 使用triangle库进行三角剖分
clc;
clear all;
load('./data/innerCircle.mat');
load('./data/outerCircle.mat');
load('./data/edgeInPlane.mat');
versInPlane = [innerCircle; outerCircle];
[~, trisInPlane] = triangle(versInPlane, edgeInPlane, mean(innerCircle), 'NoBoundarySteiners');
drawMesh(versInPlane, trisInPlane);
 
