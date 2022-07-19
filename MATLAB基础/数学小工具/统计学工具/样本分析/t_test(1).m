%单总体均值μ检验(正态总体小样本，总体标准差sigma未知)：
clc;
clear all;
alpha = 0.01;
miu0 = 1000;
%%
%(a)基本统计量直接给出：
X_bar = 98.9;
n = 16;
S = 42.3;
%%
%(b)基本统计量需要计算：
X = [774 649 1210 546 431 612];
n = numel(X);
X_bar = mean(X);
S = sqrt(var(X));

%%
t = (X_bar-miu0)/(S/sqrt(n));
disp(strcat('检测统计量t =',num2str(t)));
%%
%(a)算P值：左边检验：
p = tcdf(t,n);
disp(strcat('P值 =',num2str(p)));
%%
%(b)算P值：右边检验：
p = tcdf(t,n);
P = 1-p;
disp(strcat('P值 =',num2str(P)));
%%
%(c)算P值：双边检验：
p = tcdf(-abs(t),n);
P = 2*p;
disp(strcat('P值 =',num2str(P)));

%%
%(a)单边检验求t临界值：
beta = 1-alpha;
t_alpha = tinv(beta,n);
disp(strcat('t临界值为自由度为n的t分布上alpha分位点：',num2str(t_alpha)));
%%
%(a)双边检验求t临界值：
beta = 1-alpha/2;
t_HalfAlpha = tinv(beta,n);
disp(strcat('t临界值为自由度为n的t分布上alpha/2分位点：',num2str(t_HalfAlpha)));








%%
%匹配样本：双总体均值差miu1-miu2检验
clc;
clear all;
alpha = 0.01;
d0 = 0;                  %总体均值差的假设值。

X1 = [102 101 94 79 79];
X2 = [175 169 182 146 144];
d = X1-X2;
n = numel(X1);
X1_bar = mean(X1);
X2_bar = mean(X2);
d_bar = mean(d);
S1 = sqrt(var(X1));
S2 = sqrt(var(X2));
Sd = sqrt(var(d));

t = (d_bar-d0)/(Sd/sqrt(n));
disp(strcat('检测统计量t =',num2str(t)));

%%
%(a)算P值：左边检验：
p = tcdf(t,n-1);
disp(strcat('P值 =',num2str(p)));
%%
%(b)算P值：右边检验：
p = tcdf(t,n-1);
P = 1-p;
disp(strcat('P值 =',num2str(P)));
%%
%(c)算P值：双边检验：
p = tcdf(-abs(t),n-1);
P = 2*p;
disp(strcat('P值 =',num2str(P)));

%%
%(a)单边检验求t临界值：
beta = 1-alpha;
t_alpha = tinv(beta,n-1);
disp(strcat('自由度为(n-1)的t分布上alpha分位点为：',num2str(t_alpha)));
%%
%(a)双边检验求t临界值：
beta = 1-alpha/2;
t_HalfAlpha = tinv(beta,n-1);
disp(strcat('自由度为(n-1)的t分布上alpha/2分位点为：',num2str(t_HalfAlpha)));














%%
%正太总体小样本，sigma未知，样本容量不相等：n1≠n2,双总体均值差miu1-miu2检验
clc;
clear all;
alpha = 0.01;
d0 = 0;  

%%
%(a)基本统计量已知，不需要算
n1 = 9;
n2 = 40;
X1_bar = 70;
X2_bar = 63.2;
d_bar = X1_bar-X2_bar;
S1 = 1.5;
S2 = 2.7;
%%
%(b)基本统计量需要算：

%%
niu = round((S1^2/n1+S2^2/n2)^2/((S1^2/n1)^2/(n1-1)+(S2^2/n2)^2/(n2-1)));
t = (d_bar-d0)/sqrt(S1^2/n1+S2^2/n2);
disp(strcat('检测统计量t =',num2str(t)));

%%
%(a)算P值：左边检验：
p = tcdf(t,niu);
disp(strcat('P值 =',num2str(p)));
%%
%(b)算P值：右边检验：
p = tcdf(t,niu);
P = 1-p;
disp(strcat('P值 =',num2str(P)));
%%
%(c)算P值：双边检验：
p = tcdf(-abs(t),niu);
P = 2*p;
disp(strcat('P值 =',num2str(P)));

%%
%(a)单边检验求t临界值：
beta = 1-alpha;
t_alpha = tinv(beta,niu);
disp(strcat('自由度为niu的t分布上alpha分位点为：',num2str(t_alpha)));
%%
%(a)双边检验求t临界值：
beta = 1-alpha/2;
t_HalfAlpha = tinv(beta,niu);
disp(strcat('自由度为niu的t分布上alpha/2分位点为：',num2str(t_HalfAlpha)));














%%
%正太总体小样本，sigma未知，样本容量相等：n1=n2,双总体均值差miu1-miu2检验
clc;
clear all;
alpha = 0.05;
d0 = 0;  

%%
%(a)基本统计量已知，不需要算
n1 = 12;
n2 = 12;
n = n1;
X1_bar = -20.5;
X2_bar = -15.08333;
d_bar = X1_bar-X2_bar;
S1 = 12.38401;
S2 = 15.62317;
%%
%(b)基本统计量需要算：

%%
t = (d_bar-d0)/sqrt(S1^2/n1+S2^2/n2);
disp(strcat('检测统计量t =',num2str(t)));

%%
%(a)算P值：左边检验：
p = tcdf(t,2*n-2);
disp(strcat('P值 =',num2str(p)));
%%
%(b)算P值：右边检验：
p = tcdf(t,2*n-2);
P = 1-p;
disp(strcat('P值 =',num2str(P)));
%%
%(c)算P值：双边检验：
p = tcdf(-abs(t),2*n-2);
P = 2*p;
disp(strcat('P值 =',num2str(P)));

%%
%(a)单边检验求t临界值：
beta = 1-alpha;
t_alpha = tinv(beta,2*n-2);
disp(strcat('自由度为2*n-2的t分布上alpha分位点为：',num2str(t_alpha)));
%%
%(a)双边检验求t临界值：
beta = 1-alpha/2;
t_HalfAlpha = tinv(beta,2*n-2);
disp(strcat('自由度为2*n-2的t分布上alpha/2分位点为：',num2str(t_HalfAlpha)));









%%
%对于两个总体相关性的假设检验
%零假设：两个总体完全线性不相关，相关系数rou == 0
clc;
clear all;
%!!!数据必须是列向量
alpha = 0.05;

X1 = [102 101 94 79 79]';
X2 = [175 169 182 146 144]';
n = numel(X1);
X1_bar = mean(X1);
X2_bar = mean(X2);
S1 = sqrt(var(X1));
S2 = sqrt(var(X2));

%%
%计算样本相关系数
r = corr(X1,X2,'type','Pearson');
disp(strcat('样本相关系数r =',num2str(r)));

%%
%对总归相关性的假设检验
t = abs(r)/sqrt((1-r^2)/(n-2));      %自由度为n-2
disp(strcat('检测统计量t =',num2str(t)));

%%
%计算P值
p = tcdf(t,n-2);
P = 1-p;
disp(strcat('P值 =',num2str(P)));


