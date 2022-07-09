clc;
clear all;

%% 
%   sparse()——由稠密矩阵得到稀疏矩阵。
clc;
clear all;
m1 = [ 0   0   0   5
      0   2   0   0
      1   3   0   0
      0   0   4   0];
sm1 = sparse(m1);
rowInfo = [1,2,3];
colInfo = [3,4,5];
valueList = [0.1, 0.2, 0.3];
sm2 = sparse(rowInfo, colInfo, valueList, 5, 6);

rowInfo = [1,2,3,4,5,6];
colInfo = [1,2,3,4,5,6];
valueList = [ 1 4
              2 5
              3 6 ];          % 取元素的时候列优先
sm3 = sparse(rowInfo, colInfo, valueList, 6, 6);

rowInfo = [1,1,3,4,5,5];
colInfo = [1,1,3,4,5,5];      % 相同行列信息的值累加 
valueList = [ 1 4
              2 5
              3 6 ];        
sm4 = sparse(rowInfo, colInfo, valueList, 6, 6);

%   full()——由稀疏矩阵得到稠密矩阵。
m11 = full(sm1);
m2 = full(sm2);

% nnz()—— 返回稀疏矩阵中的非零元素数。
nonZero = nnz(sm1);
 
% nonzeros()—— 返回 稀疏矩阵的所有非零元素，存储在一个列向量内。
m3 = nonzeros(sm1);

% nzmax()—— 返回为稀疏矩阵的非零项分配的存储空间量。
size = nzmax(sm1);

% [row, col, value] = find(spM) ——获得稀疏矩阵中非零元素的信息；
elemInfo.row = [];
elemInfo.col = [];
elemInfo.value = [];
[elemInfo.row, elemInfo.col, elemInfo.value] = find(sm1);

