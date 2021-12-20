% �����������
clear;
tol = 1e-10;
f = @(x) x^3-x-1;
g = @(x) (x+1)^(1/3);
%g = @(x) x^3-1;

fprintf('��������: g=%s\n',char(g));

N = 100; % ��������������
x0 = 1.5; % ������ʼֵ
for k = 1 : N
    x = g(x0);
    fprintf('k=%2d,x=%.4e, f(x)=%.4e\n',k, x, f(x));
    if(abs(x-x0) < tol) break; end;
    x0 = x;
end
fprintf('��������: %d\n', k);