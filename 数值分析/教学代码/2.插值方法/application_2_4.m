% ���ò�ֵ��ʵ����ѧ����
clear; clc;
x = 0 : 1 : 5;
fprintf('sin����ֵ                my_sin����ֵ             ���\n');
for k = 1 : length(x)
    y = sin(x(k));
    y_star = my_sin(x(k));
    error = abs(y - y_star);    
    fprintf('%.10d       %.10d      %.10d\n', y, y_star, error);
end