%  ��ʾ���� LU �ֽ�Ľ��շ�ʽ
clear;
A=[1 2 3; 2 5 2; 3 1 5];
b=[14;18;20];

% LU���Ƿֽⷨ
[L, U, x] = my_lu(A, b)

% ����ԪLU���Ƿֽⷨ
%[L, U, x] = my_lu_with_column_pivoting(A, b)

% ֱ�ӵ���matlab��������
x_ = A\b