% ���Ϲ�ʽ����ֵ����
clear;
format long;
f = @(x) sin(x)./x; % �����������㣬ʹ�ú������������������������
a = 0; 
b = 1;
n = 8;
h = (b-a)/n;

X = a : h : b;
Y = f(X);
Y(1) = 1;  %  y(1)=0/sin(0)

% �������ι�ʽ
Tn = Y(1) + 2*sum(Y(2:n)) + Y(n+1);
Tn = Tn*h/2;
fprintf('�������ι�ʽ����Ľ���ֵΪ: %.10f\n', Tn);

% ���� Simpson ��ʽ
Sn = Y(1) + Y(n+1) + 4*sum(Y(2:2:n)) + 2*sum(Y(3:2:n-1));
Sn = Sn*h/3;
fprintf('����Simpson��ʽ����Ľ���ֵΪ: %.10f\n', Sn);

fprintf('��ȷֵΪ: %.10f\n', quad(f, 0, 1));