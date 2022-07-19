clc;
clear all;
%!!!数据必须是列向量
alpha = 0.05;

X = [102 101 94 79 79]';
Y = [175 169 182 146 144]';
n = numel(X);
X1_bar = mean(X);
X2_bar = mean(Y);
S1 = sqrt(var(X));
S2 = sqrt(var(Y));

%%
%计算样本相关系数
r = corr(X,Y,'type','Pearson');
disp(strcat('样本相关系数r =',num2str(r)));

%%
%最小二乘法拟合
%回归方程为y = beta1*x+beta0; beta1是斜率，beta0是截距。
%b1,b0分别是样本估计的斜率和截距。
b1 = (n*sum(X.*Y)-sum(X)*sum(Y))/(n*sum(X.^2)-sum(X)^2);
b0 = mean(Y)-b1*mean(X);
disp(strcat('拟合方程斜率估计值为b1 =',num2str(b1)));
disp(strcat('拟合方程截距为b0 =',num2str(b0)));

%%
%对总体相关系数ρ=0的假设检验：
%零假设：完全线性不相关，ρ=0
t = abs(r)/sqrt((1-r^2)/(n-2));      %自由度为n-2
disp(strcat('检测统计量t =',num2str(t)));

%%
%计算P值
p = tcdf(t,n-2);
P = 1-p;
disp(strcat('P值 =',num2str(P)));

%%
%t的临界值：
beta = 1-alpha/2;
t_HalfAlpha = tinv(beta,n-2);
disp(strcat('t的临界值为t分布上alpha/2分位点：',num2str(t_HalfAlpha)));

%%
%置信水平alpha下r的临界值：
r_critical = abs(t_HalfAlpha)/sqrt(t_HalfAlpha^2+n-2)





%%
%残差分析；
beta1 = 1;
beta0 = 2;
e = Y-(beta1*X+beta0);
SSE = sum(e.^2);
n = numel(e);
Se = sqrt(SSE/(n-2));





%%
%检验总体回归曲线斜率是否为0（检验x的值对于预测y值是否有意义→若斜率为0则y值是一个常数，不需要x来预测）
%零假设:斜率为0，beta1 == 0 
clc;
clear all;
t = (b1-beta1)/(Se/sqrt(Sxx));      %服从自由度为n-2的t分布



