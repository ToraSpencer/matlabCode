%%
%常数项级数求和，函数幂级数求和：
clc;
close all;
clear all;
syms k x;
ak = 3^(k+1)/(4^k);             %第k项的表达式。
s = symsum(ak,k,1,inf);
disp(s);

%%
%函数的泰勒级数展开：
clc;
close all;
clear all;
syms x;
f =power(x,-1/3);
t = taylor(f,x,8);              
disp(t);

%%
%求极限
clc;
clear all;
syms f x;
f = asin(x/(x+1));
L = limit(f,x,Inf);
disp(L);


