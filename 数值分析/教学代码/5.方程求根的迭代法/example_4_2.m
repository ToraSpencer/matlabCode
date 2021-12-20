% ��������ⷽ��
clear;
format long;
tol = 1e-5;
N = 100;
x0 = 10.5;

phi = @(x) exp(-x);
for k = 1 : N
    x1 = phi(x0);
    if abs(x1 - x0) < tol
        fprintf('���̵�����: %10.8f\n', x1);
        fprintf('��������: %d\n', k);
        break;
    end
    x0 = x1;
end
if k == N
    fprintf('��������ʧ��\n');
end