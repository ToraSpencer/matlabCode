%% 基本运算

%% 基本变换

%% 高级运算

% 范数
%       n = norm(X,p) 返回矩阵 X 的 p-范数，其中 p 为 1、2 或 Inf：
%       如果 p = 1，则 n 是矩阵的最大绝对列之和。
%       如果 p = 2，则 n 近似于 max(svd(X))。这与 norm(X) 等效。
%       如果 p = Inf，则 n 是矩阵的最大绝对行之和。
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