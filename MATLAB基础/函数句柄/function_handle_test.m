%函数句柄：

%%
%自定义函数的函数句柄：   
%               函数句柄f = @函数名；
%               函数句柄f = str2func('函数名');
clc;
close all;
clear all;

f = @findMax;
max = f(3,4);
disp('函数句柄：给函数一个别名。');
fprintf('两个数中较大的数为：%f\n\n',max)

g = str2func('findMax');
max2 = g(5,6);
fprintf('两个数中较大的数为：%f\n\n',max2)



%%
%匿名函数句柄： 函数句柄f = @(参数列表)函数表达式。
clc;
close all;
clear all;

g = @(x) x.^2+3;
disp('匿名函数句柄：快速简介地定义、使用函数。');
fprintf('g(3) = %d\n',g(3));
fprintf('g(18) = %d\n\n',g(18));

%%
%函数句柄作为另一个函数的参数：
clc;
clear all;

ezplot(@sin, [0, 2 * pi]);      %函数句柄@sin作为参数传入了函数ezplot;




%%
%速度测试：
clc;
close all;
clear all;

M = 1E6;

tic;
y = sinfun1(M);
t1 = toc;
disp('sinfun1是没有预先开辟内存的for循环');
fprintf('sinfun1()直接运行时间：t1 = %.3g\n',t1);

f1 = @() sinfun1(M);            %M是常量，这个函数里面参数列表为空。
t11 = timeit(f1);               %timeit(函数句柄)返回函数执行使用时间。
fprintf('sinfun1()使用函数句柄运行时间：t11 = %.3g\n',t11);
                                %调用函数句柄比调用函数速度快。

tic;
y = sinfun2(M);
t2 = toc;
disp('sinfun2是预先开辟内存了的for循环');              %预先开辟内存空间比临时开辟速度快。
fprintf('sinfun2()直接运行时间：t2 = %.3g\n',t2);

f2 = @() sinfun2(M);            
t21 = timeit(f2);
disp('sinfun3是将for循环向量化处理了的函数')             %向量化消除for循环加速显著。
fprintf('sinfun2()使用函数句柄运行时间：t21 = %.3g\n',t21);

tic;
y = sinfun3(M);                             
t3 = toc;
fprintf('sinfun3()直接运行时间：t3 = %.3g\n',t3);

f3 = @() sinfun3(M);            
t31 = timeit(f3);
fprintf('sinfun3()使用函数句柄运行时间：t31 = %.3g\n',t31);