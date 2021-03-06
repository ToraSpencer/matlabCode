clc;
clear all;
close all;

% 数据拟合分为两种：插值、逼近
%       插值生成的曲线或者曲面会穿过所有样本点，逼近则不一定。


%% 最小二乘拟合直线


% 1. (0, 10)×(0, 10)二维空间内的直线y = 0.5 *x + 3上取几个随机数据点，附加上随机扰动；
n = 15;                                     % 样本容量
x = 10 *rand(1,n);
y = 0.5.*x + 3; 
disturb = rand(1, n);
y = y + disturb;
plot(x, y,'b*');
hold on;

% 2. 线性方程组模型——拟合直线方程为y = k*x +b; n个样本点：(xi, yi);
%   理想情况是所有样本点都在拟合直线上——k*xi +b = yi;
figure(1)
A = [x', ones(n, 1)];
beta = y';                  % 线性方程组模型为A * X = beta, 其中X = [k; b]; 是n×2的超定线性方程组；

% 3. 求解法线方程——是n×n的恰定线性方程组：A' * A * x_bar == A' * beta 
Anew = A'*A;
betaNew = A' * beta;
x_bar = Anew \ betaNew;
k = x_bar(1);
b = x_bar(2);

% 4. 画出拟合直线
yNew = k.*x + b;
titleName = sprintf('拟合直线：y = %0.4fx + %0.4f', k, b);
title(titleName);
plot(x, yNew, 'r');

% 5. 度量余项的大小——余项欧式长度，平方误差SE，平均平方根误差RMSE
residual = b - A*x_bar;          % 余项向量；
length =  norm(residual);      % 余项向量的欧氏长度， 2范数
SE = length^2;                      % 平方误差，2范数的平方。
RMSE = SE/n;                        % 平均平方根误差
disp(['residual length == ', num2str(length)]);
disp(['SE == ', num2str(SE)]);
disp(['RMSE == ', num2str(RMSE)]);
