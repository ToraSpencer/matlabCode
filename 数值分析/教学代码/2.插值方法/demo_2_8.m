% interp1 ����
clear; clc;

xi=0:pi/5:2*pi; % ����ֵ����ֳ����ɵȾ�С����
yi=sin(xi); % ��ֵ�ڵ㴦�ĺ���ֵ

xh=0:pi/50:2*pi; % ��Ҫ��ֵ�ĵ�

% �ֶ����Բ�ֵ
subplot(2,2,1)
yh=interp1(xi,yi,xh); % ���ݲ�ֵ��������Ľ���ֵ
plot(xi,yi,'+r', xh,yh,'o-','LineWidth',1.5);
title('�ֶ����Բ�ֵ','FontSize',18)

% �ֶ���β�ֵ
subplot(2,2,2)
yh=interp1(xi,yi,xh,'nearst'); % ���ڽ���ֵ����
plot(xi,yi,'+r', xh,yh,'o-','LineWidth',1.5);
title('�ֶ���β�ֵ','FontSize',18)

% �ֶ����� Hermite
subplot(2,2,3)
yh=interp1(xi,yi,xh,'pchip'); 
plot(xi,yi,'+r', xh,yh,'o-','LineWidth',1.5);
title('�ֶ�����Hermite��ֵ','FontSize',18)

% �ֶ�����������ֵ
subplot(2,2,4)
yh=interp1(xi,yi,xh,'spline'); 
plot(xi,yi,'+r', xh,yh,'o-','LineWidth',1.5);
title('�ֶ�����������ֵ','FontSize',18)

