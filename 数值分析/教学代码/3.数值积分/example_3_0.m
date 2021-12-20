% ����˼�� �����򵥹�ʽ
clear; clc;
% �������� F
F = @(x) (+1.0146-0.88571.*(x-0.20046)+4.4375.*(x-0.20046).*(x-0.75346)-2.8282.*(x-0.20046).*(x-0.75346).*(x-1.1129)+2.8895.*(x-0.20046).*(x-0.75346).*(x-1.1129).*(x-2.3433)+0.62844*(x-0.20046).*(x-0.75346).*(x-1.1129).*(x-2.3433).*(x-1.9562)); 
%F = @(x) 1;
%F = @(x) x;
%F = @(x) x.^2;
%F = @(x) x.^3;
a = 1.0;   % ������������
b = 2.5;   % ������������
h = 0.05;  % ��ͼ����
xi = a : h : b;
yi = ones(size(xi)) .* F(xi);
xl = 0.5; xr = 3; yb= 0; yt = 5;    % ����ϵ��������

% �����ֵļ�������
subplot(2,3,1)
plot(xi,yi);
hold on;
stem(xi,yi,'k');
title('�����ֵļ�������');
axis([xl xr yb yt]); grid on;

% ����ι�ʽ
subplot(2,3,2)
plot(xi,yi);
hold on;
y = yi(1) * ones(size(xi));   % ��㺯��ֵ
stem(xi,y,'k');
title('��㷨');
axis([xl xr yb yt]); grid on;

% �Ҿ��ι�ʽ
subplot(2,3,3)
plot(xi,yi);
hold on;
y = yi(end) * ones(size(xi));  % �ҵ㺯��ֵ
stem(xi,y,'k');
title('�ҵ㷨');
axis([xl xr yb yt]); grid on;

% �о��ι�ʽ
subplot(2,3,4)
plot(xi,yi);
hold on;
y = F((a+b)/2) * ones(size(xi));  % �е㺯��ֵ
stem(xi,y,'k');
title('�е㷨');
axis([xl xr yb yt]); grid on;

% ���ι�ʽ
subplot(2,3,5)
y = yi(1)*(xi-b)/(a-b) + yi(end)*(xi-a)/(b-a);    % ���Բ�ֵ
plot(xi,yi);
hold on;
stem(xi,y,'k');
title('���ι�ʽ');
axis([xl xr yb yt]); grid on;

% �����߹�ʽ
subplot(2,3,6)
xz=(a+b)/2; yz = F(xz);   % �����������е� 
y = yi(1)*(xi-b).*(xi-xz)./((a-b)*(a-xz)) + yz*(xi-a).*(xi-b)./((xz-a).*(xz-b)) + yi(end)*(xi-a).*(xi-xz)./((b-a)*(b-xz));    % �����߲�ֵ
plot(xi,yi);
hold on;
stem(xi,y,'k');
title('�����߹�ʽ');
axis([xl xr yb yt]); grid on;
