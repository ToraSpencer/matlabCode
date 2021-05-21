%% 基本运算
% 右除/  左除\
clc;
clear all;
m1 = rand(3,3);
m2 = rand(3,3);
m3 = m1*m2;
temp1 = m3/m2;          % 右除： m1 == m3/m2 == m3 * inv(m2);
temp2 = m1\m3;          % 左除： m2 == m1\m3 == inv(m1)*m3;



%% 基本变换



%% 高级运算

%%    范数
%               n = norm(X,p) 返回矩阵 X 的 p-范数，其中 p 为 1、2 或 Inf：
%               如果 p = 1，则 n 是矩阵的最大绝对列之和。
%               如果 p = 2，则 n 近似于 max(svd(X))。这与 norm(X) 等效。
%               如果 p = Inf，则 n 是矩阵的最大绝对行之和。
v1 = round(rand(1, 3) * 10);
m1 = round(rand(2, 3) * 10);
disp(m1);
disp(v1);
norm1 = norm(m1, 1);
disp(norm1);
disp(norm(v1, 1));
norm2 = norm(m1, 2);
disp(norm2);
disp(norm(v1, 2));

%   矩阵求交：
A = [7 1 7 7 4]; 
B = [7 0 4 4 0];
[C,ia,ib] = intersect(A,B);
disp(C);            % 向量AB的共有元素
disp(ia);           % 共有元素的索引
disp(ib);           % 共有元素的索引

A = [2 2 2; 0 0 1; 1 2 3; 1 1 1; 2, 2, 2];
B = [1 2 3; 2 2 2; 2 2 0];
[C,ia,ib] = intersect(A,B,'rows');
C1 = intersect(A,B,'rows');
disp(C);            % 矩阵AB的共有行
disp(ia);           % 共有行的索引
disp(ib);           % 共有行的索引
disp(C1);


% ismember()————逐个检验m1中的每个元素在m2中是否存在。
m1 = [1,1,0,9,2,10, -1,1,15,1,3,1];
m2 = [1,2,3; 4, 5, 6; 7, 8, 9];
disp(ismember(m1, m2));

% union()——矩阵的并，会自动去除重复的行
m1 = [1,2,3; 4,5,6; 1,2,3];
m2 = [-1,-2,-3; 4, 5, 6; 7, 8, 9; -4, -5, -6];
m3 = union(m1, m2, 'rows');
disp(m3);


% setdiff()——矩阵元素的差集，会自动去除重复的行
m1 = [1,2,3; 4,5,6; 1,2,3];
m2 = [-1,-2,-3; 4, 5, 6; 7, 8, 9; -4, -5, -6];
m3 = setdiff(m1, m2, 'rows');
disp(m3);



%% lu分解——[L,U,P,Q,D] = lu(S)是高斯消元法的矩阵形式
% S == D * inv(P) * L * U * inv(Q)
% L是下三角矩阵，U是上三角矩阵，inv(P)是左乘置换矩阵，inv(Q)是右乘置换矩阵，D是对角缩放矩阵。
clc;
clear all;
m1 = [1,2,3,4;5,6,7,8;9,10,0,12;13,14,15,0];
disp(rank(m1));
disp(det(m1));
sm1 = sparse(m1);
[sL, sU, sP, sQ, sD] = lu(sm1);
L = full(sL);
U = full(sU);
P = full(sP);
Q = full(sQ);
D = full(sD);
disp(m1);
disp(D * inv(P) * L * U * inv(Q) );
disp(D);



%% 稀疏矩阵

%   sparse()——由稠密矩阵得到稀疏矩阵。
%%
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

%%


%   full()——由稀疏矩阵得到稠密矩阵。
m11 = full(sm1);
m2 = full(sm2);

% nnz()—— 返回稀疏矩阵中的非零元素数。
nonZero = nnz(sm1);
 
% nonzeros()—— 返回 稀疏矩阵的所有非零元素，存储在一个列向量内。
m3 = nonzeros(sm1);

% nzmax()—— 返回为稀疏矩阵的非零项分配的存储空间量。
size = nzmax(sm1);



%%
clc;
clear all;
m1 = [11,12,13,14;21,22,23,24; 31,32,33,34];
v1 = [1,2,3];
v2 = [3,4];
m2 = m1(v1, v2);
disp(m2);

V = v1'*v2;
disp(V);
