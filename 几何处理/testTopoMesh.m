clc;
clear all;
close all;


%% 顶点邻接矩阵
clc;
clear all;
close all;
[vers1, tris1,Q] = cube(2, 2, 2);
[vers2, tris2] = readOBJ('./data/holeTooth.obj');
writeOBJ('holeTooth.obj', vers2, tris2);
versCount1 = size(vers1, 1);
versCount2 = size(vers2, 1);

edges1 = [tris1(:,2) tris1(:,3); tris1(:,3) tris1(:,1); tris1(:,1) tris1(:,2)];    % 有向边，与三角片正法向对应；
tmpSM = sparse(edges1(:,1), edges1(:,2), 1);   
tmpM = full(tmpSM);
adjM1 = tmpM > 0;                % 顶点有向邻接矩阵(logical)；下标(i, j)的元素为true == 存在有向边(i, j);

edges2 = [tris2(:,2) tris2(:,3); tris2(:,3) tris2(:,1); tris2(:,1) tris2(:,2)];
tmpSM = sparse(edges2(:,1), edges2(:,2), 1);   
adjSM2 = tmpSM > 0; 

% 无孤立点无重复三角片的三角网格是封闭 == 顶点有向邻接矩阵是一个对称矩阵；
nonDlM1 = adjM1 - adjM1';               % 非双向(double linked)边矩阵；若下标为(i, j)的元素不为0，则表示(i, j)为非双向边；
nonDlSM2 = adjSM2 - adjSM2';       
flag1 = all(all(nonDlM1 == 0));
flag2 = all(all(nonDlSM2 == 0));
disp(strcat('adjM1是否是对称矩阵？网格1是否封闭？——', num2str(flag1)));
disp(strcat('adjSM2是否是对称矩阵？网格2是否封闭？——', num2str(flag2)));


%% 通过顶点邻接矩阵找洞：
% 通过有向邻接矩阵计算非封闭网格的边缘——边缘处的边，其正向和逆向边并不成对存在；
elemInfo.row = [];
elemInfo.col = [];
elemInfo.value = [];
[elemInfo.row, elemInfo.col, elemInfo.value] = find(nonDlSM2);
bdryEdges = [elemInfo.row, elemInfo.col];
objWriteEdges('bdryEdges.obj', bdryEdges, vers2);

% 
elemInfo22.row = [];
elemInfo22.col = [];
elemInfo22.value = [];
adjSM22 = 2* adjSM2 - adjSM2';          % 1表示是双向边，0表示无边，2和-1表示单向边；
[elemInfo22.row, elemInfo22.col, elemInfo22.value] = find(adjSM22);
disp('finished.');
 



 
