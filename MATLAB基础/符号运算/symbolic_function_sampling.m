%符号函数采样：符号函数可以视为连续函数，要获得离散的函数值数据，需要进行采样。

%%
%
clc; close all;clear all;
syms t f1 f2 f3                               %声明符号变量；

f1 = dirac(t);                                %dirac函数δ(t)=dirac(t)
f2 = sinc(t);
f3 = sin(2*pi*t).*heaviside(t);               %阶跃函数step(t)=heaviside(t)

%   选取区间观察函数取值: subs()
