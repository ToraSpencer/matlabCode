%%
%���������ͣ������ݼ�����ͣ�
clc;
close all;
clear all;
syms k x;
ak = 3^(k+1)/(4^k);             %��k��ı��ʽ��
s = symsum(ak,k,1,inf);
disp(s);

%%
%������̩�ռ���չ����
clc;
close all;
clear all;
syms x;
f =power(x,-1/3);
t = taylor(f,x,8);              
disp(t);

%%
%����
clc;
clear all;
syms f x;
f = asin(x/(x+1));
L = limit(f,x,Inf);
disp(L);


