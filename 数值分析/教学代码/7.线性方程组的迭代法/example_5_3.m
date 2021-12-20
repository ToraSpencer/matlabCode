% ���ɳڷ���������ʽ�� ��ⷽ����
clear;
% ����ֵ
A = [10, -1, -2; -1, 10, -2; -1, -1, 5];
b = [7.2; 8.3; 4.2];
tol = 1e-5;
N = 100;   % ����������
omega = 1.55;  % �ɳ�����,ȡֵ����(0, 2)
x = [0; 0; 0];
x_backup = [0; 0; 0];   
y = [0; 0; 0];
%
A_ = A;
for i = 1 : length(A)
    A_(i,i) = 0;
end

for k = 0 : N
    for j = 1 : length(A)
        y(j) = (b(j) - A_(j,:) * x) / A(j,j);
        y(j) = (1 - omega) * x(j) + omega * y(j);
        x_backup(j) = x(j);   % ���ݡ���ֵ��
        x(j) = y(j);     % ����ֵ�� �滻 ����ֵ��
    end
    if (max(abs(x_backup - y)) < tol)
        fprintf('��������: %d\n', k);
        fprintf('������ĸ�: %10.8f\n', y);
        break;
    end
end
if k == N
    fprintf('��������ʧ��\n');
end