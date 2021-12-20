% Newton �������ľֲ�������
clear; clc
f = @(x) x^3 - x -1;
df = @(x) 3*x^2 - 1;

N = 20;  % ����������
tol = 1e-6;
%x0 = 1.5; % ��ֵ�ӽ���  f(x0)*f(x)�Ķ��׵��� > 0 , ��&��ֵ���� or ͹&��ֵΪ����
x0 = 0.008; % ��ֵԶ���
for k = 1 : N
    x = x0 - f(x0)/df(x0);
    fprintf('k=%d, x=%.8f\n',k,x);
    if abs(x-x0)<tol
        fprintf('��������: %d\n', k);
        fprintf('���̵�����: %10.8f\n', x);
        break; 
    end
    x0 = x;
end
if k == N
    fprintf('��������ʧ��\n');
end
