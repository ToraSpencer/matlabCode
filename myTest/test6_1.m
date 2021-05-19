clc
close all;
clear all;

addpath('mesh process');
addpath('mixFE');

%% 6.1 计算变形所需的参数——biharm_factor_system()
load('N0.mat');
load('N1.mat');
load('omega.mat');
load('mergedToothVers.mat');
load('newTris.mat');

interest = [omega N0 N1 ];

temp = newTris( ...
ismember(newTris(:,1),interest) & ...
ismember(newTris(:,2),interest) & ...
ismember(newTris(:,3),interest),:);     % 三个顶点都在intere里的所有三角片。

temp2 = temp';

all = [N0 omega];
S = cotmatrix(mergedToothVers, newTris);  
elems = nonzeros(S);

% 求协方差矩阵
i1 = temp2(1,:); i2 = temp2(2,:); i3 = temp2(3,:); 
vers1 = mergedToothVers(i3,:) - mergedToothVers(i2,:);  
vers2 = mergedToothVers(i1,:) - mergedToothVers(i3,:); 
vers3 = mergedToothVers(i2,:) - mergedToothVers(i1,:);
vers4  = cross(vers1,vers2, 2);

temp = (vers4').^2;
temp = sum(temp);
temp = sqrt(temp);
dblA = temp';


cot12 = -dot(vers1,vers2,2)./dblA/2; 
cot23 = -dot(vers2,vers3,2)./dblA/2; 
cot31 = -dot(vers3,vers1,2)./dblA/2;
diag1 = -cot12-cot31; 
diag2 = -cot12-cot23; 
diag3 = -cot31-cot23;
rowInfo = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
colInfo = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
valueList = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];
S2 = sparse(rowInfo, colInfo, valueList, size(mergedToothVers,1), size(mergedToothVers,1));
elems2 = nonzeros(S2); 

disp('finished.');

 
 
