function [ y ] = ill_si( t, x, lambda)
%ILL_SIR Summary of this function goes here
%   ������� x ����һ��������Ϊ��Infective �ı���
% lambda   �սӴ���

y = [lambda * x(1) * (1 - x(1))]';

end

