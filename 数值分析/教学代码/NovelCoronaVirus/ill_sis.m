% ����SISģ��΢�ַ���
function [ y ] = ill_sis( t, x, lambda, miu)
%ILL_SIS Summary of this function goes here
%   ������� x ����һ��������Ϊ��Infective �ı���
% lambda   �սӴ��ʣ���ÿ������ÿ���Ⱦ������
%  miu     �������ʣ���ÿ������������ռ���������ı���

y = [lambda * x(1) * (1 - x(1)) - miu * x(1)]';

end

