%% 快捷键
% 折叠代码    CTRL + =   
% 取消折叠    CTRL + SHIFT + =
% 多行注释    CTRL + R
% 多行反注释  CTRL + T
% 打开变量、函数 CTRL + D
%           打开函数指的是打开函数文件。
%           变量只能在debug的时候打开，可以显示数据的详情。

%% 排序算法——sort()
clc;
clear all;
x = round(100 * rand(1, 10));
disp(x);
[u, v] = sort(x);
disp(u);
disp(v);

 