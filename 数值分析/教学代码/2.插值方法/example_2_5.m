% p23 ��5 ���ؽ��𲽲�ֵ
clear;
format long;
x = 0.3 :0.1 : 0.7;
y = [0.29850, 0.39646, 0.49311, 0.58813, 0.68122];
xx = 0.462; % ��ֵ��
fij(:,1) = y;
for i = 1 : 3
    for j = i+1 : 5
       fij(j,i+1) = fij(i,i)*(xx-x(j))/(x(i)-x(j)) + fij(j,i)*(xx-x(i))/(x(j)-x(i))  % ��ʾ�м��������ֵ���
    end
end
fij