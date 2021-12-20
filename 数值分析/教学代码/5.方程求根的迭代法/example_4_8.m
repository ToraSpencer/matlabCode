% �����ҽط���ⷽ��
clear;
format long;
tol = 1e-5;
N = 100;
x0 = 0.5;
x1 = 0.6;
%f = @(x) x * exp(x) - 1; % f(x)���ʽ
f = @(x) x^3 - x - 1;  %f(x)���ʽ

for k = 1 : N
    x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
    if abs(x2 - x1) < tol
        fprintf('��������: %d\n', k);
        fprintf('���̵�����: %10.8f\n', x1);
        break;
    end
    x0 = x1;
    x1 = x2;
end
if k == N
    fprintf('��������ʧ��\n');
end