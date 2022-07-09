%% 基本运算
% 右除/  左除\
clc;
clear all;
m1 = rand(3,3);
m2 = rand(3,3);
m3 = m1*m2;
temp1 = m3/m2;          % 右除： m1 == m3/m2 == m3 * inv(m2);
temp2 = m1\m3;          % 左除： m2 == m1\m3 == inv(m1)*m3;


%% 矩阵数据的增删查改：

% min(), max()——查找最值；
clc;
clear all;
vec1 = [8.1, 7.3, 6.2, 5.5, 4.4, 3.1, 2.5, 1.1];
[value, index] = min(vec1);
 
% [row, col, value] = find(M) ——获得矩阵中非零元素的信息；
M1 =  [0, 8.1, 7.3; 6.2, 0, 4.4; -3.1, 2.5, 1.1];
elemInfo.row = [];
elemInfo.col = [];
elemInfo.value = [];
[elemInfo.row, elemInfo.col, elemInfo.value] = find(M1);

% [row, col, value] = find(元素条件) —— 获得矩阵中满足条件的元素的信息；
% [row, col] = find(元素条件)  
[elemInfo.row, elemInfo.col, elemInfo.value] = find(M1 < 0);        % 这里的value是逻辑值1


%% 使用左除解线性方程组：
% 3x+y = 1; 2x-y = 2; 
clc;
clear all;
A = [3,1; 2,-1];
b = [1;2];              % A*x = b;
x1 = inv(A)*b;
x2 = A\b;               % A矩阵可以是稀疏矩阵。
disp(x1);
disp(x2);



%% 基本变换



%% 高级运算



%%    范数
%               n = norm(X), 返回矩阵X的2范数，即sqrt(x1^2+x2^2+x3^2 ....);
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

%   矩阵求交：intersect()——
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



%% 奇异值分解——svd()
clc;
clear all;

%       [U,S,V] = svd(A) 执行矩阵 A 的奇异值分解，因此 A = U*S*V'，其中S是对角矩阵，元素为A的奇异值。
A = [1 0 1; -1 -2 0; 0 1 -1];
[U,S,V] = svd(A); 
disp(U);
disp(S);
disp(V);

%       s = svd(A) 以降序顺序返回矩阵 A 的奇异值。


%% 索引矩阵：
clc;
clear all;

% 矩阵的条件表达式返回一个索引矩阵，其中满足条件的元素为1，不满足的为0；
M1 =  [0, 8.1, 7.3; 6.2, 0, 4.4; -3.1, 2.5, 1.1];
IM1 = (M1 < 0);     
M2 = M1.*IM1;               % 与逻辑矩阵点乘，可以筛选元素；

 
