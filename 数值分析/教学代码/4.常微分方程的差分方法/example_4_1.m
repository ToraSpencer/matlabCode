% ŷ��������ֵ����
clear;
f = @(x,y) (y - 2 * x / y);
h = 0.2;   % ����
x = 0:h:2;  % ��ɢ��
y(1) = 1;   % ��ֵ
y_xn = (1 + 2 * x).^(1/2);   % ������
plot(x,y_xn,'k-'); 
hold on;
% ŷ������
for i = 2 : length(x)
    y(i) = y(i-1) + h * f(x(i-1),y(i-1));
    plot(x(i-1:i), y(i-1:i), 'o-');
    hold on;
end
% plot(x,y,'o-',x,y_xn,'k-'); 
legend('��ȷ��', 'ŷ����');
