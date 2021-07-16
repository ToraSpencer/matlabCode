clc;
clear all;
close all;
%%
[vers, tris] = readOBJ('./data/mountain.obj');

% 画出边界线
plot_boundary_orange(vers, tris);


% 寻找边缘三角片
flag = on_boundary(tris);       % flag向量
bdryTriIdx = find(flag);            % flag向量转化为索引向量

% 计算边缘曲线长度
length = boundary_length(vers, tris);