%%
%用符号描写函数；符号描述可以直接写出表达式，然后任意地截取区间来观察函数。
clc 
clear all
syms t x y z1 z2 z3 z4                            %声明符号变量；
x=sin(2*pi*t).*heaviside(t);                %符号函数；
y=exp(-t);                                  %阶跃函数step(t)=heaviside(t)
z1=x+y;
z2=x*y;
z3=dirac(t);                                %dirac函数δ(t)=dirac(t)
z4=sinc(t);
%%
%选取区间观察函数取值：两种方法：使用subs或eval命令。
x=subs(x,t,[-1:0.05:2]);                   %用[-1:0.05:2]替换掉函数x中的t，得到x函数值的分布，但是x的数据类型仍然为符号;
y=subs(y,t,[-1:0.05:2]);
t=[-1:0.05:2];
z1=eval(z1);                               %eval()作用是将括号内的字符串视为语句并运行。得到的z1数据类型是double;
z2=eval(z2);
z3=eval(z3);
%%
%画图：两种方法：(1)符号转为double→生成图框→plot命令;(2)直接对符号和符号函数t,z4使用ezplot命令。
x=double(x);                              %把符号x转换为double数值；
y=double(y); 
z1=double(z1); 
z2=double(z2); 

figure;                                   %生成一个图框；
hold on;
subplot(1,2,1);                           %图框分成一行两列的两块，对编号为1的那块作图；
plot(t,x,'r',t,y,'b-h');                  %x图线选为红色，y图线选为蓝色六角形线；
xlabel('t(s)') ;                          %填写x轴说明。
ylabel('signals');
title('x,y');                             %填写标题
legend('x','y')                           %填写图例  
subplot(1,2,2);
plot(t,z1,'g-o',t,z2,t,z3,'y');
xlabel('t(s)')  ;                         
ylabel('signals');
title('z1,z2,z3');
legend('z1','z2','z3') ;

figure;                                   %再生成一个新的图框
plot(t,x,'k',t,y,'b',t,z1,'g-o',t,z2,'r',t,z3,'y');
legend('x','y','z1','z2','z3');   

figure;
ezplot(z4,[-10,10]); 




