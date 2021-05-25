clc;
close all;
clear all;
%%
% bi_L, bi_U, bi_P, bi_Q, bi_D, bi_S, finalVers, Omega, N0, N1)
load('bi_L.mat');
load('bi_U.mat');
load('bi_P.mat');
load('bi_Q.mat');
load('bi_D.mat');
load('bi_S.mat');
load('A.mat');
load('mergedToothVers.mat');
load('mergeRegionVers.mat');
load('patientTransform.mat');
load('omega.mat');
load('N0.mat');
load('N1.mat');
load('mergeRegionIdx.mat');

finalVers = mergedToothVers;
versCount = length(mergedToothVers);
ptfCount = length(patientTransform);
mrVersCount = length(mergeRegionVers);

%% 找出patientTransform中距离mergeRegion中某一点最近的点 
for i = 1:mrVersCount
   ver = mergeRegionVers(i,:);
   temp = repmat(ver, ptfCount, 1) - patientTransform;
   temp = temp.^2;
   temp = sum(temp, 2); 
   Distance = sqrt(temp);       % mergeRegion中当前考察的点到patientTransform点的距离。
   
   [u,v]=sort(Distance);
    minDisIndex = v(1);        % distance中最小元素的索引, 考虑用KD树来做？？？？
  
   finalVers(mergeRegionIdx(i),:) = patientTransform(minDisIndex,:);
end

    all = [N0 omega];
    allSize = size(all,2);
  
%% 解稀疏矩阵A的线性方程组Ax = b;
for i = 1:3

    % b向量
    temp1 = -bi_S(all,N0) * finalVers(N0,i);
    temp2 = bi_S(all,N1) * finalVers(N1,i);
  rhs_Dx = temp1 - temp2;
  rhs_Dy = zeros(size(omega,2),1);
  rhs = [ rhs_Dx; rhs_Dy];   

    % 原程序使用QPLU矩阵来计算稀疏矩阵表示的线性方程组，matlab只能解稠密矩阵表示的线性方程组。
  Afull = full(A);
  sol = linsolve(Afull, rhs);
  temp = sol(allSize+1:end); 
  finalVers(omega,i) = sol(allSize+1:end); 
end

OBJwriteVertices('最终顶点.obj', finalVers);

disp('finished.');
