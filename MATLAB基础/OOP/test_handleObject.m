clc;
clear all;
handle1 = vector2D(2,7);
handle2 = handle1;
%%
%继承handle的类实例化出来的对象本质上是句柄，两个句柄指向了同一个value对象，调用一个句柄的方法改变了这个value对象的属性，另一个句柄对象的属性也随之变化。
handle1.x = 5;
disp(handle1.x);
disp(handle2.x);

%%
%handle对象调用析构函数：
handle1.delete();
fprintf('%d   %d\n',isvalid(handle1),isvalid(handle2));  %handle所指的数据被释放，所有handle对象无效化。

%%
%无效化后的handle对象可以重新指向数据。
handle1 = vector2D(2,7);
handle2 = handle1;
disp(handle1);          
disp(handle2);

%%
%对handle对象使用clear函数。
clear handle1;
fprintf('%d\n',isvalid(handle2)); 
