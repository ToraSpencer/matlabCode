% ��������ⷽ�̵�Ψһ����
clear;
format long;
tol = 1e-6;
N = 20;
x0 = 80.5;

for k = 1 : N
    x1 = (x0 + 1) ^ (1 / 3);
    if abs(x1 - x0) < tol
        fprintf('���̵�����: %10.8f\n', x1);
        break;
    end
    x0 = x1;
end
if k == N
    fprintf('��������ʧ��');
end