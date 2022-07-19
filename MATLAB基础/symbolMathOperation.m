%%
%�÷�����д������������������ֱ��д�����ʽ��Ȼ������ؽ�ȡ�������۲캯����
clc 
clear all
syms t x y z1 z2 z3 z4                            %�������ű�����
x=sin(2*pi*t).*heaviside(t);                %���ź�����
y=exp(-t);                                  %��Ծ����step(t)=heaviside(t)
z1=x+y;
z2=x*y;
z3=dirac(t);                                %dirac������(t)=dirac(t)
z4=sinc(t);
%%
%ѡȡ����۲캯��ȡֵ�����ַ�����ʹ��subs��eval���
x=subs(x,t,[-1:0.05:2]);                   %��[-1:0.05:2]�滻������x�е�t���õ�x����ֵ�ķֲ�������x������������ȻΪ����;
y=subs(y,t,[-1:0.05:2]);
t=[-1:0.05:2];
z1=eval(z1);                               %eval()�����ǽ������ڵ��ַ�����Ϊ��䲢���С��õ���z1����������double;
z2=eval(z2);
z3=eval(z3);
%%
%��ͼ�����ַ�����(1)����תΪdouble������ͼ���plot����;(2)ֱ�ӶԷ��źͷ��ź���t,z4ʹ��ezplot���
x=double(x);                              %�ѷ���xת��Ϊdouble��ֵ��
y=double(y); 
z1=double(z1); 
z2=double(z2); 

figure;                                   %����һ��ͼ��
hold on;
subplot(1,2,1);                           %ͼ��ֳ�һ�����е����飬�Ա��Ϊ1���ǿ���ͼ��
plot(t,x,'r',t,y,'b-h');                  %xͼ��ѡΪ��ɫ��yͼ��ѡΪ��ɫ�������ߣ�
xlabel('t(s)') ;                          %��дx��˵����
ylabel('signals');
title('x,y');                             %��д����
legend('x','y')                           %��дͼ��  
subplot(1,2,2);
plot(t,z1,'g-o',t,z2,t,z3,'y');
xlabel('t(s)')  ;                         
ylabel('signals');
title('z1,z2,z3');
legend('z1','z2','z3') ;

figure;                                   %������һ���µ�ͼ��
plot(t,x,'k',t,y,'b',t,z1,'g-o',t,z2,'r',t,z3,'y');
legend('x','y','z1','z2','z3');   

figure;
ezplot(z4,[-10,10]); 




