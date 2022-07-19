%计算特殊的几种离散型随机变量某点的概率，累计概率，数字特征。
%%
%二项分布：
clc;
clear all;
p = 0.7;
n = 15;
k = 6;
P1 = binopdf(k,n,p);
P2 = binocdf(k,n,p);
disp(sprintf('概率为%g的基本事件在%d重伯努利实验中发生%d次的概率为%g',p,n,k,P1));
disp(sprintf('概率为%g的基本事件在%d重伯努利实验中至少发生%d次的概率为%g',p,n,k,P2));
Ex = n*p;
Dx = n*p*(1-p);
disp(sprintf('期望为%g，方差为%g',Ex,Dx));


