%%
%һԪ�����ֵļ��㣺
clc;
clear all;
syms f x 
f = x^2;
a = 1;
b = 2;
S = int(f,x,a,b);           %���غ���f(x)������[a,b]�ϵĻ��֡�
disp(S);

%%
%һԪ�������ֵļ��㣺
clc;
clear all;
syms f x 
f = x^2;
S = int(f,x);               %���غ���f(x)�Ĳ������֡�
disp(S);

%%
%һԪ����һ��΢�ֵļ��㣺
clc;
clear all;
syms f x 
f = -5*exp(-x);
D = diff(f,'x');
disp(D);

%%
%һԪ����n��΢�ֵļ��㣺
clc;
clear all;
syms f x 
n = 9;
f = exp(-x);
D = diff(f,'x',n);
disp(D);