% ���ٵ���������ⷽ�� ȱ������Ҫ��⵼��ֵL
clear;
format long;
tol = 1e-5;
N = 100;
x0 = 0.5;
phi = @(x) exp(-x);
L = -exp(-x0);  % ��ֵ��ĵ���ֵ

for k = 1 : N
    x1 = (1/(1-L)) * (phi(x0) - L * x0);
    if abs(x1 - x0) < tol
        fprintf('���̵�����: %10.8f\n', x1);
        break;
    end
    x0 = x1;
end
if k == N
    fprintf('��������ʧ��\n');
end