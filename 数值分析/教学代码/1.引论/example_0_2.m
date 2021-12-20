% ���ַ��󷽳̽��Ƹ�
clear;
f = @(x) x ^ 3 - x -1;
a = 1;
b = 1.5;
tol = 1e-15;
n = 1;

y0 = f(a);
while (abs(b - a) > tol)
    x = (a + b) / 2;
    y = f(x);
    if (y * y0 > 0)
        a = x;
    else
        b = x;
    end
    n = n + 1;
end
x;
fprintf(' x = %.10d\n', x);
fprintf('��������: %d\n', n);