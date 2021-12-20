% Newton��ɽ�� ����ֵѡ�񲻵�ʱ ������ɽ����lamda
clear;
format long;
tol = 1e-5;
N = 100;
x0 = 0.6;
lamda = 1;
f = @(x) x^3 - x - 1;  %f(x)���ʽ
df = @(x) 3*x^2 - 1;
fprintf('f(x)�ĳ�ֵ: %d\n', abs(f(x0)));
for k = 1 : N
    x1 = x0 - f(x0)/ df(x0);
     while abs(f(x1)) > abs(f(x0))
         lamda = lamda / 2;   % ��ɽ���Ӽ���
         fprintf('��ɽ���Ӽ�����ֵ: %d\n', lamda);
         x1 = x0 - lamda * f(x0) / df(x0);
     end
    fprintf('���ε���f(x)��ֵ: %d\n', abs(f(x1)));
    if abs(x1 - x0) < tol
        fprintf('��������: %d\n', k);
        fprintf('���̵�����: %10.8f\n', x1);
        break;
    end
    x0 = x1;
end
if k == N
    fprintf('��������ʧ��\n');
end