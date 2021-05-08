%% 基本运算

%% 基本变换

%% 高级运算

%% 
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