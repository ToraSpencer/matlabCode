%单总体均值μ检验：
clc;
clear all;

alpha = 0.01;
miu0 = 1000;
sigma = 0.62;
%%
%(a)估计值直接给出：
X_bar = 98.2;
n = 106;
S = 1;
%%
%(b)估计值需要计算：
X = [774 649 1210 546 431 612];
n = numel(X);
X_bar = mean(X);
S = sqrt(var(X));

%%
%(a)总体标准差sigma已知：
Z = (X_bar-miu0)/(sigma/sqrt(n));
disp(strcat('检测统计量Z =',num2str(Z)));
%%
%(b)大样本情形下，总体标准差sigma未知：
Z = (X_bar-miu0)/(S/sqrt(n));
disp(strcat('检测统计量Z =',num2str(Z)));

%%
%(a)算P值：左边检验：
p = normcdf(Z);
disp(strcat('P值 =',num2str(p)));
%%
%(b)算P值：右边检验：
p = normcdf(Z);
P = 1-p;
disp(strcat('P值 =',num2str(P)));
%%
%(c)算P值：双边检验：
p = normcdf(-abs(Z));
P = 2*p;
disp(strcat('P值 =',num2str(P)));

%%
beta = 1-alpha;
Z_alpha = norminv(beta);
disp(strcat('标准正态分布上alpha分位点为：',num2str(Z_alpha)));












%%
%单总体比例p检验：
clc;
clear all;
alpha = 0.01;

%%
%(a)估计值直接给出：
p_hat = 0.2;
n = 1501;
p0 = 0.25;
%%
%(b)估计值需要计算：
n = 514;
n1 = 236;
p_hat = n1/n;
p0 = 0.5;

%%
Z = (p_hat-p0)/sqrt(p0*(1-p0)/n);
disp(strcat('检测统计量Z =',num2str(Z)));

%%
%(a)算P值：左边检验：
p = normcdf(Z);
disp(strcat('P值 =',num2str(p)));
%%
%(b)算P值：右边检验：
p = normcdf(Z);
P = 1-p;
disp(strcat('P值 =',num2str(P)));
%%
%(c)算P值：双边检验：
p = normcdf(-abs(Z));
P = 2*p;
disp(strcat('P值 =',num2str(P)));

%%
%(a)单边检验求Z临界值：
beta = 1-alpha;
Z_alpha = norminv(beta);
disp(strcat('Z临界值为标准正态分布上alpha分位点：',num2str(Z_alpha)));
%%
%(a)双边检验求Z临界值：
beta = 1-alpha/2;
Z_HalfAlpha = norminv(beta);
disp(strcat('Z临界值为标准正态分布上alpha/2分位点：',num2str(Z_HalfAlpha)));









%%
%双总体比例差p1-p2检验
clc;
clear all;
alpha = 0.01;
d0 = 0;  

%%
%(a)基本统计量已知，不需要算
%%
%(b)基本统计量需要算：
n1 = 15;
N1 = 343;
n2 = 27;
N2 = 294;
p1_hat = n1/N1;
p2_hat = n2/N2;

%%
Z = (p1_hat-p2_hat-d0)/sqrt(p1_hat*(1-p1_hat)/N1+p2_hat*(1-p2_hat)/N2);
disp(strcat('检测统计量Z =',num2str(Z)));

%%
%(a)算P值：左边检验：
p = normcdf(Z);
disp(strcat('P值 =',num2str(p)));
%%
%(b)算P值：右边检验：
p = normcdf(Z);
P = 1-p;
disp(strcat('P值 =',num2str(P)));
%%
%(c)算P值：双边检验：
p = normcdf(-abs(Z));
P = 2*p;
disp(strcat('P值 =',num2str(P)));

%%
%(a)单边检验求Z临界值：
beta = 1-alpha;
Z_alpha = norminv(beta);
disp(strcat('标准正态分布上alpha分位点为：',num2str(Z_alpha)));
%%
%(a)双边检验求Z临界值：
beta = 1-alpha/2;
Z_HalfAlpha = norminv(beta);
disp(strcat('标准正态分布上alpha/2分位点为：：',num2str(Z_HalfAlpha )));

