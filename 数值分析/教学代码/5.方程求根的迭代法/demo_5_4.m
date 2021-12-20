% ����Newton���Ŀ�����ʽ
clear; clc

c = 115;
N = 10; % ����������
tol = 1e-6;
x0 = 15; % ��ֵ
for k = 1 : N
    x = 0.5 * (x0 + c / x0);
    fprintf('k=%d, x=%.8f\n',k,x);
    if abs(x - x0) < tol
        fprintf('���̵�����: %10.8f\n', x);
        break;
    end
    x0 = x;
end
fprintf('��������: %d\n', k);
