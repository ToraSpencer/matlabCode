%���ź������������ź���������Ϊ����������Ҫ�����ɢ�ĺ���ֵ���ݣ���Ҫ���в�����

%%
%
clc; close all;clear all;
syms t f1 f2 f3                               %�������ű�����

f1 = dirac(t);                                %dirac������(t)=dirac(t)
f2 = sinc(t);
f3 = sin(2*pi*t).*heaviside(t);               %��Ծ����step(t)=heaviside(t)

%   ѡȡ����۲캯��ȡֵ: subs()
