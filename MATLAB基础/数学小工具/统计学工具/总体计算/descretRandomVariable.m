%��������ļ�����ɢ���������ĳ��ĸ��ʣ��ۼƸ��ʣ�����������
%%
%����ֲ���
clc;
clear all;
p = 0.7;
n = 15;
k = 6;
P1 = binopdf(k,n,p);
P2 = binocdf(k,n,p);
disp(sprintf('����Ϊ%g�Ļ����¼���%d�ز�Ŭ��ʵ���з���%d�εĸ���Ϊ%g',p,n,k,P1));
disp(sprintf('����Ϊ%g�Ļ����¼���%d�ز�Ŭ��ʵ�������ٷ���%d�εĸ���Ϊ%g',p,n,k,P2));
Ex = n*p;
Dx = n*p*(1-p);
disp(sprintf('����Ϊ%g������Ϊ%g',Ex,Dx));


