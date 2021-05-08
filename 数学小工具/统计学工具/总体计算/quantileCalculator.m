clc;
close all;
clear all;
%%
%���׼��̬�ֲ�����alpha��λ��
alpha = 0.01;
beta = 1-alpha;
z = norminv(beta);    %��׼��̬�ۻ��ֲ������ķ�����
disp('��׼��̬�ֲ�����alpha��Ϊ��Ϊ��')
disp(z);

%%
%�󿨷��ֲ��ϵ�alpha��λ�㣺
alpha = 0.05;
beta = 1-alpha;
n = 15;
chi_alpha = chi2inv(beta,n);
disp('���ɶ�Ϊn�Ŀ����ֲ���alpha��λ��Ϊ��');
disp(chi_alpha);

%%
%��t�ֲ���alpha��λ�㣺
alpha = 0.01;
beta = 1-alpha;
n = 15;
t_alpha = tinv(beta,n-1);
disp('���ɶ�Ϊn��t�ֲ���alpha��λ��Ϊ��');
disp(t_alpha);

%%
%��F�ֲ���alpha��λ��
alpha = 0.05;
beta = 1-alpha;
n1 = 9;
n2 = 1;
F_alpha = finv(beta,n1,n2);
disp('���ɶ�Ϊ(n1,n2)��F�ֲ���alpha��λ��Ϊ��');
disp(F_alpha);

