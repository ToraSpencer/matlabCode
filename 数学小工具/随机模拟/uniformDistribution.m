clc;
clear all;
p1 = unidpdf(2,6);
cp1 = unidcdf(2,6);
disp(strcat('掷骰子实验中，出现2的概率为',num2str(p1)));
disp(strcat('掷骰子实验中，2的累计概率为',num2str(cp1)));

M = unidrnd(6,[1,10]);
disp('掷骰子掷10次，结果为：');
M


data = unidrnd(6,[1,100000]);
hist(data,[1,2,3,4,5,6]);           %掷骰子100000次实验结果的频数分布直方图。