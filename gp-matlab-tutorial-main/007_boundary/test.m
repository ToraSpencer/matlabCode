clc;
clear all;
close all;
%%
[vers, tris] = readOBJ('./data/mountain.obj');

% �����߽���
plot_boundary_orange(vers, tris);


% Ѱ�ұ�Ե����Ƭ
flag = on_boundary(tris);       % flag����
bdryTriIdx = find(flag);            % flag����ת��Ϊ��������

% �����Ե���߳���
length = boundary_length(vers, tris);