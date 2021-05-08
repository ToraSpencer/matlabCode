%用样本统计量估计总体分布参数miu,sigma的范围区间。
%%
%1. 估计miu,sigma已知
clc;
clear all;
sigma = 495;
alpha = 0.05;
%%
%(a)基本统计量已知，不需要计算。
n = 75;
X_bar = 3433;
%%
%(b)基本统计量需要计算。
X = [62 61 61 57 61 54 59 58 59 69 60 67];
X_bar = mean(X);
n = numel(X);
%%
halfAlpha = alpha/2;
beta = 1-halfAlpha;
z_halfAlpha = norminv(beta);
error = sigma/sqrt(n)*z_halfAlpha;
leftMargin = X_bar-error;
rightMargin = X_bar+error;
disp(strcat('误差为：',num2str(error)));
disp(strcat('置信区间为：','（',num2str(leftMargin),',',num2str(rightMargin),')'));








%%
%2. 估计miu,sigma未知
clc;
clear all;
alpha = 0.01; 
%%
%(a)基本统计量已知，不需要计算。
n = 20;
X_bar = 17.38;
S = 8.24;
%%
%(b)基本统计量需要计算。
X = [62 61 61 57 61 54 59 58 59 69 60 67];
S = sqrt(var(X));
X_bar = mean(X);
n = numel(X);
%%
halfAlpha = alpha/2;
beta = 1-halfAlpha;
t_halfAlpha = tinv(beta,n-1);
error = S/sqrt(n)*t_halfAlpha;
leftMargin = X_bar-error;
rightMargin = X_bar+error;
disp(strcat('误差为：',num2str(error)));
disp(strcat('置信区间为：','（',num2str(leftMargin),',',num2str(rightMargin),')'));












%%
%3. 估计sigma
clc;
clear all;
alpha = 0.05;
halfAlpha = alpha/2;
%%
%(a)基本统计量已知，不需要计算。
S = 1;
n = 10;
%%
%(b)基本统计量需要计算。
X = [62 61 61 57 61 54 59 58 59 69 60 67];
S = sqrt(var(X));
n = numel(X);
%%
leftMargin = sqrt((n-1)*S^2/chi2inv(1-halfAlpha,n-1));
rightMargin = sqrt((n-1)*S^2/chi2inv(halfAlpha,n-1));
disp(strcat('置信区间为：','（',num2str(leftMargin),',',num2str(rightMargin),')'));

