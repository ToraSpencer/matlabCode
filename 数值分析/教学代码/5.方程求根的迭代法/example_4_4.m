% ���ؽ𷽷����ٵ���������ⷽ��
clear;
format long;
tol = 1e-10;
N = 100;
x0 = 0.5;
phi = @(x) exp(-x);
%phi = @(x) (x+1)^(1/3);

for k = 1 : N
    x1 = phi(x0);
    x2 = phi(x1);
    x2 = x2 - (x2 - x1)^2 / (x2 - 2 * x1 + x0);
    if abs(x2 - x0) < tol
        fprintf('���̵�����: %10.8f\n', x1);
        break;
    end
    x0 = x2;
end
if k == N
    fprintf('��������ʧ��\n');
end
fprintf('��������: %d\n', k);