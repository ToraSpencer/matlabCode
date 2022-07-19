%约束条件分三种：
%1.决策变量的边界 lb(1)<x<ub(1);  lb(2)<y<ub(2)
%2.决策变量的等式组 
%3.决策变量的不等式组。
%linprog函数只能求目标函数的最小值。要求最大值的问题可以转化为求z = -f'*x的最小值。

%x:线性规划中的决策变量，如[x y]'，[x1,x2,x3,x4....]'
%z = f'*x:目标函数。f是目标函数的系数行向量。
%A*x<=b：约束条件中的线性不等式组。其中A为系数矩阵，b为右边值组成的列向量。
%Aeq*x = beq:约束条件中的等式组。
%lb：约束条件中的决策变量的下界组成的列向量。
%ub：约束条件中的决策变量的上界组成的列向量。
clc;
clear all;
%%
%x = f(f,A,b)。当没有最小值的时候，返回的x为空矩阵。
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
    disp('目标函数没有最小值。');
end
if numel(x_maxValueProblem)==0
    disp('目标函数没有最大值。');
end
    
if numel(x_minValueProblem)~=0
    z = f*x_minValueProblem;
    disp(strcat('目标函数最小值为',num2str(f*x_minValueProblem)));
    disp('对应的决策变量坐标点为：');
    disp(x_minValueProblem');
end

if numel(x_maxValueProblem)~=0
    z = f*x_maxValueProblem;
    disp(strcat('目标函数最大值为',num2str(f*x_maxValueProblem)));
    disp('对应的决策变量坐标点为：');
    disp(x_maxValueProblem');
end



