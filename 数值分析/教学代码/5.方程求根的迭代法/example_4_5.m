% Newton���������߷�����ⷽ��
clear;
format long;
tol = 1e-5;
N = 100;
x0 = 0.5;
phi = @(x) x - (x - exp(-x)) / (1 + x);  % ��ɷ��̶�Ӧ��ţ�ٹ�ʽ

for k = 1 : N
    x1 = phi(x0);
    if abs(x1 - x0) < tol
        fprintf('���̵�����: %10.8f\n', x1);
        break;
    end
    x0 = x1;
end
if k == N
    fprintf('��������ʧ��\n');
end