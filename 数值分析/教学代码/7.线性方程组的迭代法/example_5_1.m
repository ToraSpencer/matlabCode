% Jacobi ��������������ʽ����ⷽ����
clear;
% ����ֵ
A = [10, -1, -2; -1, 10, -2; -1, -1, 5];
b = [7.2; 8.3; 4.2];
tol = 1e-5;
x = [0; 0; 0];
y = [0; 0; 0];

A_ = A;
for i = 1 : length(A)
    A_(i,i) = 0;    % �Խ���Ԫ������Ϊ0
end
for i = 0 : 19
%     fprintf('��%d�ε���: \n', i);
    
%     y = (b - A_ * x) ./ diag(A);  % Jacobi������ʽ(������ʽ)

    for j = 1 : length(A)
        y(j,1) = (b(j) - sum(A_(j,:)*x))/A(j,j);  % Jacobi������ʽ(������ʽ)
    end

    if (max(abs(x - y)) < tol)
        fprintf('��������: %d\n', i);
        fprintf('������ĸ�: %10.8f\n', y);
        break;
    end
    x = y;
end