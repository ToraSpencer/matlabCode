%计算特殊的几种连续型随机变量的概率密度，累计概率，数字特征。

%%
%正态分布
clc;
clear all;
mu = 63.6;
sigma = 2.5;
x1 = 66.5;
x2 = 71.5;
y1 = normcdf(x1,mu,sigma);
y2 = normcdf(x2,mu,sigma);
disp('X大于x1的概率为');
disp(1-y1);
disp('X小于x1的概率为')
disp(y1);
disp('X大于x2的概率为')
disp(1-y2);
disp('X小于x2的概率为')
disp(y2);
disp('X在x1和x2之间的概率为');
disp(y2-y1)

%%
%卡方分布：
clc;
clear all;
n = 3;
x = 4.069;
p = chi2cdf(x,n);
disp('X<x的概率为：');
disp(p);

%%
%t分布：
clc;
clear all;
n = 23;
x1 = 2.618;
x2 = 2.618;
y1 = tcdf(x1,n);
y2 = tcdf(x2,n);
disp('X大于x1的概率为')
disp(1-y1);
disp('X小于x1的概率为')
disp(y1);
disp('X大于x2的概率为')
disp(1-y2);
disp('X小于x2的概率为')
disp(y2);
disp('X在x1和x2之间的概率为');
disp(y2-y1)

%%
%F分布：
clc;
clear all;
n1 = 9;
n2 = 6;
x1 = 2.618;
x2 = 2.618;
y1 = fcdf(x1,n1,n2);
y2 = fcdf(x2,n1,n2);
disp('X大于x1的概率为')
disp(1-y1);
disp('X小于x1的概率为')
disp(y1);
disp('X大于x2的概率为')
disp(1-y2);
disp('X小于x2的概率为')
disp(y2);
disp('X在x1和x2之间的概率为');
disp(y2-y1)
