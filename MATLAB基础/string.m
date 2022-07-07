clc;
clear all;
close all;


%% 字符串拼接，格式化字符串. strcat(), num2str();
str1 = 'aaabbccc';
f1 = 3.1415;
str2 = strcat('str1', num2str(f1));
disp(str2);
