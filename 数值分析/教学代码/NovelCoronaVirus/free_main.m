% ��������ģ����ʾ����

ts = 0 : 20;   % ʱ������
lambda = 0.3;  % ÿ������ÿ���Ⱦ����
x0 = 1;   % ��ʼ������
infective = x0 * exp(lambda * ts);   % ������ָ������
plot(ts, infective);
title('��������ģ��');
grid;