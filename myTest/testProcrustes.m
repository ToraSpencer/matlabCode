clc;
clear all;
load('movedRootTooth1.mat');
load('point11.mat');
load('point2.mat');


%% 1.
n  = size(point11,1);
 
% Center at the origin.
center1 = mean(point11,1);
center2 = mean(point2,1);
X0 = point11 - repmat(center1, n, 1);
Y0 = point2 - repmat(center2, n, 1);
 
normX = sqrt(trace(X0*X0'));
normY = sqrt(trace(Y0*Y0'));
X0 = X0 / normX;
Y0 = Y0 / normY;
 
A = X0' * Y0;
[U, ~, V] = svd(A);
T = V * U';
 
tc = center1 - center2*T;
tcc = repmat(tc, size(movedRootTooth1, 1), 1);
movedRootTooth22 = movedRootTooth1 + tcc;
OBJwriteVertices('movedRootTooth22.obj', movedRootTooth22);
disp('finished.');