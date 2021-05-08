%�������ֵ�̼��飺
clc;
clear all;

alpha = 0.01;
miu0 = 1000;
sigma = 0.62;
%%
%(a)����ֱֵ�Ӹ�����
X_bar = 98.2;
n = 106;
S = 1;
%%
%(b)����ֵ��Ҫ���㣺
X = [774 649 1210 546 431 612];
n = numel(X);
X_bar = mean(X);
S = sqrt(var(X));

%%
%(a)�����׼��sigma��֪��
Z = (X_bar-miu0)/(sigma/sqrt(n));
disp(strcat('���ͳ����Z =',num2str(Z)));
%%
%(b)�����������£������׼��sigmaδ֪��
Z = (X_bar-miu0)/(S/sqrt(n));
disp(strcat('���ͳ����Z =',num2str(Z)));

%%
%(a)��Pֵ����߼��飺
p = normcdf(Z);
disp(strcat('Pֵ =',num2str(p)));
%%
%(b)��Pֵ���ұ߼��飺
p = normcdf(Z);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));
%%
%(c)��Pֵ��˫�߼��飺
p = normcdf(-abs(Z));
P = 2*p;
disp(strcat('Pֵ =',num2str(P)));

%%
beta = 1-alpha;
Z_alpha = norminv(beta);
disp(strcat('��׼��̬�ֲ���alpha��λ��Ϊ��',num2str(Z_alpha)));












%%
%���������p���飺
clc;
clear all;
alpha = 0.01;

%%
%(a)����ֱֵ�Ӹ�����
p_hat = 0.2;
n = 1501;
p0 = 0.25;
%%
%(b)����ֵ��Ҫ���㣺
n = 514;
n1 = 236;
p_hat = n1/n;
p0 = 0.5;

%%
Z = (p_hat-p0)/sqrt(p0*(1-p0)/n);
disp(strcat('���ͳ����Z =',num2str(Z)));

%%
%(a)��Pֵ����߼��飺
p = normcdf(Z);
disp(strcat('Pֵ =',num2str(p)));
%%
%(b)��Pֵ���ұ߼��飺
p = normcdf(Z);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));
%%
%(c)��Pֵ��˫�߼��飺
p = normcdf(-abs(Z));
P = 2*p;
disp(strcat('Pֵ =',num2str(P)));

%%
%(a)���߼�����Z�ٽ�ֵ��
beta = 1-alpha;
Z_alpha = norminv(beta);
disp(strcat('Z�ٽ�ֵΪ��׼��̬�ֲ���alpha��λ�㣺',num2str(Z_alpha)));
%%
%(a)˫�߼�����Z�ٽ�ֵ��
beta = 1-alpha/2;
Z_HalfAlpha = norminv(beta);
disp(strcat('Z�ٽ�ֵΪ��׼��̬�ֲ���alpha/2��λ�㣺',num2str(Z_HalfAlpha)));









%%
%˫���������p1-p2����
clc;
clear all;
alpha = 0.01;
d0 = 0;  

%%
%(a)����ͳ������֪������Ҫ��
%%
%(b)����ͳ������Ҫ�㣺
n1 = 15;
N1 = 343;
n2 = 27;
N2 = 294;
p1_hat = n1/N1;
p2_hat = n2/N2;

%%
Z = (p1_hat-p2_hat-d0)/sqrt(p1_hat*(1-p1_hat)/N1+p2_hat*(1-p2_hat)/N2);
disp(strcat('���ͳ����Z =',num2str(Z)));

%%
%(a)��Pֵ����߼��飺
p = normcdf(Z);
disp(strcat('Pֵ =',num2str(p)));
%%
%(b)��Pֵ���ұ߼��飺
p = normcdf(Z);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));
%%
%(c)��Pֵ��˫�߼��飺
p = normcdf(-abs(Z));
P = 2*p;
disp(strcat('Pֵ =',num2str(P)));

%%
%(a)���߼�����Z�ٽ�ֵ��
beta = 1-alpha;
Z_alpha = norminv(beta);
disp(strcat('��׼��̬�ֲ���alpha��λ��Ϊ��',num2str(Z_alpha)));
%%
%(a)˫�߼�����Z�ٽ�ֵ��
beta = 1-alpha/2;
Z_HalfAlpha = norminv(beta);
disp(strcat('��׼��̬�ֲ���alpha/2��λ��Ϊ����',num2str(Z_HalfAlpha )));

