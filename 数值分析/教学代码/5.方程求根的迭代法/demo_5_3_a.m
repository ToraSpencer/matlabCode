% Newton �������ĳ�ֵѡ��
clear; clc
f = @(x) x.^3 - x -1;

%x0 = 1.5; % ��ֵ�ӽ���  f(x0)*f(x)�Ķ��׵��� > 0 , ��&��ֵ���� or ͹&��ֵΪ����
%x0 = 0.5; % ��ֵԶ���

x=linspace(-2,2,300);
y=f(x);

[y0,i]=sort(abs(y));    % ����i�������Ԫ�ص����
plot(x,y,x(i(1)),y0(1),'rp');
ha=gca;
set(ha,'ygrid','on')

