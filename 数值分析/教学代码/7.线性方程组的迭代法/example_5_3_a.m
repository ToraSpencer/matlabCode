% ���ɳڷ���������ʽ�� ��ⷽ����
clear;
% ����ֵ
A = [10, -1, -2; -1, 10, -2; -1, -1, 5];
b = [7.2; 8.3; 4.2];
tol = 1e-5;
N = 100;   % ����������
omega = 1.07;  % �ɳ�����,ȡֵ����(0, 2),ȡֵ1.0xʱ����Ч����
x = [0; 0; 0];

D = diag(diag(A));  % A �ĶԽ��߲���
L = D - tril(A);    % -L Ϊ A ���ϸ������ǲ���
U = D - triu(A);    % -U Ϊ A ���ϸ������ǲ���

for k = 0 : N
    y = (D-omega*L) \ ( ((1-omega)*D + omega*U)*x + omega*b );
    if (max(abs(x - y)) < tol)
        fprintf('��������: %d\n', k);
        fprintf('������ĸ�: %10.8f\n', y);
        break;
    end
    x = y;
end
if k == N
    fprintf('��������ʧ��\n');
end