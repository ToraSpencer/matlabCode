%%
%一元定积分的计算：
clc;
clear all;
syms f x 
f = x^2;
a = 1;
b = 2;
S = int(f,x,a,b);           %返回函数f(x)在区间[a,b]上的积分。
disp(S);

%%
%一元不定积分的计算：
clc;
clear all;
syms f x 
f = x^2;
S = int(f,x);               %返回函数f(x)的不定积分。
disp(S);

%%
%一元函数一阶微分的计算：
clc;
clear all;
syms f x 
f = -5*exp(-x);
D = diff(f,'x');
disp(D);

%%
%一元函数n阶微分的计算：
clc;
clear all;
syms f x 
n = 9;
f = exp(-x);
D = diff(f,'x',n);
disp(D);