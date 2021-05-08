clc;
close all;
clear all;
%%
%求标准正态分布的上alpha分位点
alpha = 0.01;
beta = 1-alpha;
z = norminv(beta);    %标准正态累积分布函数的反函数
disp('标准正态分布的上alpha分为点为：')
disp(z);

%%
%求卡方分布上的alpha分位点：
alpha = 0.05;
beta = 1-alpha;
n = 15;
chi_alpha = chi2inv(beta,n);
disp('自由度为n的卡方分布上alpha分位点为：');
disp(chi_alpha);

%%
%求t分布上alpha分位点：
alpha = 0.01;
beta = 1-alpha;
n = 15;
t_alpha = tinv(beta,n-1);
disp('自由度为n的t分布上alpha分位点为：');
disp(t_alpha);

%%
%求F分布上alpha分位点
alpha = 0.05;
beta = 1-alpha;
n1 = 9;
n2 = 1;
F_alpha = finv(beta,n1,n2);
disp('自由度为(n1,n2)的F分布上alpha分位点为：');
disp(F_alpha);

