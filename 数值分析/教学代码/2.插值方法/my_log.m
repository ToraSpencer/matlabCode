function [ y ] = my_log( xl )
%MY_LOG �Բ�ֵ����ʵ����Ȼ��������
%   Detailed explanation goes here
b = 1 : 4;
yb = log(b);
[m, c] = d_d(b, yb);
y = nest(c, xl, b);
end

