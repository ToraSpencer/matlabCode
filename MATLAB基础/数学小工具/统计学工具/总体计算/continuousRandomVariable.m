%��������ļ�����������������ĸ����ܶȣ��ۼƸ��ʣ�����������

%%
%��̬�ֲ�
clc;
clear all;
mu = 63.6;
sigma = 2.5;
x1 = 66.5;
x2 = 71.5;
y1 = normcdf(x1,mu,sigma);
y2 = normcdf(x2,mu,sigma);
disp('X����x1�ĸ���Ϊ');
disp(1-y1);
disp('XС��x1�ĸ���Ϊ')
disp(y1);
disp('X����x2�ĸ���Ϊ')
disp(1-y2);
disp('XС��x2�ĸ���Ϊ')
disp(y2);
disp('X��x1��x2֮��ĸ���Ϊ');
disp(y2-y1)

%%
%�����ֲ���
clc;
clear all;
n = 3;
x = 4.069;
p = chi2cdf(x,n);
disp('X<x�ĸ���Ϊ��');
disp(p);

%%
%t�ֲ���
clc;
clear all;
n = 23;
x1 = 2.618;
x2 = 2.618;
y1 = tcdf(x1,n);
y2 = tcdf(x2,n);
disp('X����x1�ĸ���Ϊ')
disp(1-y1);
disp('XС��x1�ĸ���Ϊ')
disp(y1);
disp('X����x2�ĸ���Ϊ')
disp(1-y2);
disp('XС��x2�ĸ���Ϊ')
disp(y2);
disp('X��x1��x2֮��ĸ���Ϊ');
disp(y2-y1)

%%
%F�ֲ���
clc;
clear all;
n1 = 9;
n2 = 6;
x1 = 2.618;
x2 = 2.618;
y1 = fcdf(x1,n1,n2);
y2 = fcdf(x2,n1,n2);
disp('X����x1�ĸ���Ϊ')
disp(1-y1);
disp('XС��x1�ĸ���Ϊ')
disp(y1);
disp('X����x2�ĸ���Ϊ')
disp(1-y2);
disp('XС��x2�ĸ���Ϊ')
disp(y2);
disp('X��x1��x2֮��ĸ���Ϊ');
disp(y2-y1)
