%  ��ʾ���� ��ȥ��
clear;
A=[1,1,1; 0,4,-1; 2,-2,1];
b=[6;5;1];

%  ���� ��˹��ȥ��
%x = my_ge(A, b)

%  ���� ����Ԫ��˹��ȥ��
x = my_ge_with_column_pivoting(A, b)

% ֱ�ӵ���matlab��������
x_ = A\b

