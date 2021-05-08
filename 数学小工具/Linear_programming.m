%Լ�����������֣�
%1.���߱����ı߽� lb(1)<x<ub(1);  lb(2)<y<ub(2)
%2.���߱����ĵ�ʽ�� 
%3.���߱����Ĳ���ʽ�顣
%linprog����ֻ����Ŀ�꺯������Сֵ��Ҫ�����ֵ���������ת��Ϊ��z = -f'*x����Сֵ��

%x:���Թ滮�еľ��߱�������[x y]'��[x1,x2,x3,x4....]'
%z = f'*x:Ŀ�꺯����f��Ŀ�꺯����ϵ����������
%A*x<=b��Լ�������е����Բ���ʽ�顣����AΪϵ������bΪ�ұ�ֵ��ɵ���������
%Aeq*x = beq:Լ�������еĵ�ʽ�顣
%lb��Լ�������еľ��߱������½���ɵ���������
%ub��Լ�������еľ��߱������Ͻ���ɵ���������
clc;
clear all;
%%
%x = f(f,A,b)����û����Сֵ��ʱ�򣬷��ص�xΪ�վ���
f = [20 60];
A = [-90 -60;-2 -2;-5 -10];
b = [-360 -10 -30];
b = b';
Aeq = [];
beq = [];
beq = beq';
lb = [0 0]; 
ub = [Inf Inf]; 
x_minValueProblem = linprog(f,A,b,Aeq,beq,lb,ub);
x_maxValueProblem = linprog(-f,A,b,Aeq,beq,lb,ub);

if numel(x_minValueProblem)==0
    disp('Ŀ�꺯��û����Сֵ��');
end
if numel(x_maxValueProblem)==0
    disp('Ŀ�꺯��û�����ֵ��');
end
    
if numel(x_minValueProblem)~=0
    z = f*x_minValueProblem;
    disp(strcat('Ŀ�꺯����СֵΪ',num2str(f*x_minValueProblem)));
    disp('��Ӧ�ľ��߱��������Ϊ��');
    disp(x_minValueProblem');
end

if numel(x_maxValueProblem)~=0
    z = f*x_maxValueProblem;
    disp(strcat('Ŀ�꺯�����ֵΪ',num2str(f*x_maxValueProblem)));
    disp('��Ӧ�ľ��߱��������Ϊ��');
    disp(x_maxValueProblem');
end



