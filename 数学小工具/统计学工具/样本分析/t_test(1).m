%�������ֵ�̼���(��̬����С�����������׼��sigmaδ֪)��
clc;
clear all;
alpha = 0.01;
miu0 = 1000;
%%
%(a)����ͳ����ֱ�Ӹ�����
X_bar = 98.9;
n = 16;
S = 42.3;
%%
%(b)����ͳ������Ҫ���㣺
X = [774 649 1210 546 431 612];
n = numel(X);
X_bar = mean(X);
S = sqrt(var(X));

%%
t = (X_bar-miu0)/(S/sqrt(n));
disp(strcat('���ͳ����t =',num2str(t)));
%%
%(a)��Pֵ����߼��飺
p = tcdf(t,n);
disp(strcat('Pֵ =',num2str(p)));
%%
%(b)��Pֵ���ұ߼��飺
p = tcdf(t,n);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));
%%
%(c)��Pֵ��˫�߼��飺
p = tcdf(-abs(t),n);
P = 2*p;
disp(strcat('Pֵ =',num2str(P)));

%%
%(a)���߼�����t�ٽ�ֵ��
beta = 1-alpha;
t_alpha = tinv(beta,n);
disp(strcat('t�ٽ�ֵΪ���ɶ�Ϊn��t�ֲ���alpha��λ�㣺',num2str(t_alpha)));
%%
%(a)˫�߼�����t�ٽ�ֵ��
beta = 1-alpha/2;
t_HalfAlpha = tinv(beta,n);
disp(strcat('t�ٽ�ֵΪ���ɶ�Ϊn��t�ֲ���alpha/2��λ�㣺',num2str(t_HalfAlpha)));








%%
%ƥ��������˫�����ֵ��miu1-miu2����
clc;
clear all;
alpha = 0.01;
d0 = 0;                  %�����ֵ��ļ���ֵ��

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
disp(strcat('���ͳ����t =',num2str(t)));

%%
%(a)��Pֵ����߼��飺
p = tcdf(t,n-1);
disp(strcat('Pֵ =',num2str(p)));
%%
%(b)��Pֵ���ұ߼��飺
p = tcdf(t,n-1);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));
%%
%(c)��Pֵ��˫�߼��飺
p = tcdf(-abs(t),n-1);
P = 2*p;
disp(strcat('Pֵ =',num2str(P)));

%%
%(a)���߼�����t�ٽ�ֵ��
beta = 1-alpha;
t_alpha = tinv(beta,n-1);
disp(strcat('���ɶ�Ϊ(n-1)��t�ֲ���alpha��λ��Ϊ��',num2str(t_alpha)));
%%
%(a)˫�߼�����t�ٽ�ֵ��
beta = 1-alpha/2;
t_HalfAlpha = tinv(beta,n-1);
disp(strcat('���ɶ�Ϊ(n-1)��t�ֲ���alpha/2��λ��Ϊ��',num2str(t_HalfAlpha)));














%%
%��̫����С������sigmaδ֪��������������ȣ�n1��n2,˫�����ֵ��miu1-miu2����
clc;
clear all;
alpha = 0.01;
d0 = 0;  

%%
%(a)����ͳ������֪������Ҫ��
n1 = 9;
n2 = 40;
X1_bar = 70;
X2_bar = 63.2;
d_bar = X1_bar-X2_bar;
S1 = 1.5;
S2 = 2.7;
%%
%(b)����ͳ������Ҫ�㣺

%%
niu = round((S1^2/n1+S2^2/n2)^2/((S1^2/n1)^2/(n1-1)+(S2^2/n2)^2/(n2-1)));
t = (d_bar-d0)/sqrt(S1^2/n1+S2^2/n2);
disp(strcat('���ͳ����t =',num2str(t)));

%%
%(a)��Pֵ����߼��飺
p = tcdf(t,niu);
disp(strcat('Pֵ =',num2str(p)));
%%
%(b)��Pֵ���ұ߼��飺
p = tcdf(t,niu);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));
%%
%(c)��Pֵ��˫�߼��飺
p = tcdf(-abs(t),niu);
P = 2*p;
disp(strcat('Pֵ =',num2str(P)));

%%
%(a)���߼�����t�ٽ�ֵ��
beta = 1-alpha;
t_alpha = tinv(beta,niu);
disp(strcat('���ɶ�Ϊniu��t�ֲ���alpha��λ��Ϊ��',num2str(t_alpha)));
%%
%(a)˫�߼�����t�ٽ�ֵ��
beta = 1-alpha/2;
t_HalfAlpha = tinv(beta,niu);
disp(strcat('���ɶ�Ϊniu��t�ֲ���alpha/2��λ��Ϊ��',num2str(t_HalfAlpha)));














%%
%��̫����С������sigmaδ֪������������ȣ�n1=n2,˫�����ֵ��miu1-miu2����
clc;
clear all;
alpha = 0.05;
d0 = 0;  

%%
%(a)����ͳ������֪������Ҫ��
n1 = 12;
n2 = 12;
n = n1;
X1_bar = -20.5;
X2_bar = -15.08333;
d_bar = X1_bar-X2_bar;
S1 = 12.38401;
S2 = 15.62317;
%%
%(b)����ͳ������Ҫ�㣺

%%
t = (d_bar-d0)/sqrt(S1^2/n1+S2^2/n2);
disp(strcat('���ͳ����t =',num2str(t)));

%%
%(a)��Pֵ����߼��飺
p = tcdf(t,2*n-2);
disp(strcat('Pֵ =',num2str(p)));
%%
%(b)��Pֵ���ұ߼��飺
p = tcdf(t,2*n-2);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));
%%
%(c)��Pֵ��˫�߼��飺
p = tcdf(-abs(t),2*n-2);
P = 2*p;
disp(strcat('Pֵ =',num2str(P)));

%%
%(a)���߼�����t�ٽ�ֵ��
beta = 1-alpha;
t_alpha = tinv(beta,2*n-2);
disp(strcat('���ɶ�Ϊ2*n-2��t�ֲ���alpha��λ��Ϊ��',num2str(t_alpha)));
%%
%(a)˫�߼�����t�ٽ�ֵ��
beta = 1-alpha/2;
t_HalfAlpha = tinv(beta,2*n-2);
disp(strcat('���ɶ�Ϊ2*n-2��t�ֲ���alpha/2��λ��Ϊ��',num2str(t_HalfAlpha)));









%%
%����������������Եļ������
%����裺����������ȫ���Բ���أ����ϵ��rou == 0
clc;
clear all;
%!!!���ݱ�����������
alpha = 0.05;

X1 = [102 101 94 79 79]';
X2 = [175 169 182 146 144]';
n = numel(X1);
X1_bar = mean(X1);
X2_bar = mean(X2);
S1 = sqrt(var(X1));
S2 = sqrt(var(X2));

%%
%�����������ϵ��
r = corr(X1,X2,'type','Pearson');
disp(strcat('�������ϵ��r =',num2str(r)));

%%
%���ܹ�����Եļ������
t = abs(r)/sqrt((1-r^2)/(n-2));      %���ɶ�Ϊn-2
disp(strcat('���ͳ����t =',num2str(t)));

%%
%����Pֵ
p = tcdf(t,n-2);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));


