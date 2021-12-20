% Newton���������߷������ÿ�����ʽ��115��ƽ����
clear;
format long;
tol = 1e-6;
N = 100;
x0 = 10;
c = 115;
phi = @(x) (x + c / x) / 2;  % ���η��̶�Ӧ��ţ�ٹ�ʽ

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