clc;
clear all;
%%
%Ԥ�ȿ����ڴ�ռ����ʱ�����ٶȿ졣
%���ú�������ȵ��ú����ٶȿ졣
%����������forѭ������������
%%
M = 1E6;

tic;
y = sinfun1(M);
t1 = toc;
fprintf('t1 = %.3g\n',t1);


f1 = @() sinfun1(M);            %�޲��������������
t11 = timeit(f1);               %timeit(�������)���غ���ִ��ʹ��ʱ�䡣
fprintf('t11 = %.3g\n',t11);

tic;
y = sinfun2(M);
t2 = toc;
fprintf('t2 = %.3g\n',t2);

f2 = @() sinfun2(M);            
t21 = timeit(f2);
fprintf('t21 = %.3g\n',t21);

tic;
y = sinfun3(M);
t3 = toc;
fprintf('t3 = %.3g\n',t3);

f3 = @() sinfun3(M);            
t31 = timeit(f3);
fprintf('t31 = %.3g\n',t31);




