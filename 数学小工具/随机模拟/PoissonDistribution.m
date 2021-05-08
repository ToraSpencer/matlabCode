clc;
clear all;
x= 0:15;
p = poisspdf(x,5);
cp = poisscdf(x,5);

subplot(2,1,1);
bar(x,p);
xlim([-1,16]);

subplot(2,1,2);
stairs(x,cp);
%%
data = poissrnd(5,[1,10000]);
x1 = 0:1:30;
hist(data,x1);