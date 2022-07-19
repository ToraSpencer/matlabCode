clc;
clear all;
%%
%预先开辟内存空间比临时开辟速度快。
%调用函数句柄比调用函数速度快。
%向量化消除for循环加速显著。
%%
M = 1E6;

tic;
y = sinfun1(M);
t1 = toc;
fprintf('t1 = %.3g\n',t1);


f1 = @() sinfun1(M);            %无参匿名函数句柄。
t11 = timeit(f1);               %timeit(函数句柄)返回函数执行使用时间。
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




