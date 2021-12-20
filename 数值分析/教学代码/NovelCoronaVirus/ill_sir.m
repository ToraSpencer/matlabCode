% SIRģ��΢�ַ���
function [ y ] = ill_sir( t, x )
%ILL_SIR Summary of this function goes here
%   ������� x ������������������Ϊ��Suspectable��Infective �ı���

lambda = 0.2586;  % �սӴ���
miu = 0.0821;  % �������ʣ���ÿ������������ռ���������ı���

y=[lambda * x(1) * x(2) - miu * x(1), -lambda * x(1) * x(2)]';

end

