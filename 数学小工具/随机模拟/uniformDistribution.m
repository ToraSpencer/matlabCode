clc;
clear all;
p1 = unidpdf(2,6);
cp1 = unidcdf(2,6);
disp(strcat('������ʵ���У�����2�ĸ���Ϊ',num2str(p1)));
disp(strcat('������ʵ���У�2���ۼƸ���Ϊ',num2str(cp1)));

M = unidrnd(6,[1,10]);
disp('��������10�Σ����Ϊ��');
M


data = unidrnd(6,[1,100000]);
hist(data,[1,2,3,4,5,6]);           %������100000��ʵ������Ƶ���ֲ�ֱ��ͼ��