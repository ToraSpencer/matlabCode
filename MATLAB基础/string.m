clc;
clear all;
close all;


%% �ַ���ƴ�ӣ���ʽ���ַ���. strcat(), num2str();
str1 = 'aaabbccc';
f1 = 3.1415;
str2 = strcat('str1', num2str(f1));
disp(str2);
