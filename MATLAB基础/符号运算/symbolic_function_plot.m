%符号函数画图

%%
%
clc; close all;clear all;
syms t f
f = sin(2*pi*t)*t;                %符号函数；

%   画图：直接对符号t和符号函数f(t)使用ezplot()
figure; ezplot(f,[-10,10]); title('f(t) = sin(2*pi*t)*t');




