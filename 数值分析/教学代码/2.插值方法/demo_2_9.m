% spline ����
clear; clc;

xi=[27.7, 28, 29, 30]; % ��ֵ�ڵ�
yi=[4.1, 4.3, 4.1, 3.0]; % �ڵ㴦�ĺ���ֵ
df0=3.0; dfn=-4.0; % �߽�����
pp=spline(xi,[df0, yi, dfn]);

xh=27.7:0.1:30; % ��Ҫ��ֵ�ĵ�
yh=ppval(pp,xh); % ͨ����ֵ��õĽ���ֵ 
plot(xi,yi,'r+',xh,yh,'o-b','LineWidth',1.5,'MarkerSize',12);
