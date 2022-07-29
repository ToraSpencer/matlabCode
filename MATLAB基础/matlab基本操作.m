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


%% 画图
clc;
clear all;
close all;
x = 1:100;
y = 0.1*(x-35).^2 + 10;

figure('Name', 'xy graph');     

% plot()——画折线图
plot(x, y, 'b-*');
xlabel('axis X');
ylabel('axis Y');
hold on                  % 保留当前坐标区中的绘图，从而使新添加到坐标区中的绘图不会删除现有绘图;

% scatter() ——画散点图
x1 = 10*[1: 10];
y1 = -0.1*(x1-40).^2 + 30;
y2 = y1-30;
scatter(x1, y1, 'r');                       % 画的点默认是空心圆；
scatter(x1, y2, 'g', 'filled');             % filled指定画的点为实心圆；

% surf()——画曲面：
% meshgrid()——生成栅格
% xlabel(), ylabel——设定轴信息
[X,Y]=meshgrid(-2:0.1:2, -3: 0.1 :3);
Z=exp(-X.^2-Y.^2);
figure(2);
xlabel('x');
ylabel('y');
%           修改图像句柄的属性：
surfHandle = surf(X,Y,Z);          
surfHandle.EdgeColor = 'none';              % 隐藏边；
        


%% 读写.dat文件
clc;
clear all;
load('./data/elliSampleVers.mat');
save ./data/elliSampleVers2.mat sampleVers;
x = sampleVers(:, 1);
y = sampleVers(:, 2);
figure('Name', 'xy graph');
plot(x, y, 'b-*');


%% 调用控制台程序：
clc;
clear all;
fileName = 'G:/gitRepositories/matlabCode/MATLAB基础/data/tooth.obj';
tic
command = ['E:/workstation/smarteeproj/Exe/Win32_Release/planeCutMesh.exe', ...
   ' planeCutMesh ', fileName, ' [0, 0, 0] [0, 0, 1]'];
[status, result] = system( command );
fprintf('planeCutMesh calculation takes %f s time.\n', toc);


 