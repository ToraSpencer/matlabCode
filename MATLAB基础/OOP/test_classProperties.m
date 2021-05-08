clc;
clear all;
obj = N_randomNumber();

%%
a = obj.cumulativeProbability;      %属性没有设默认值，构造方法中又没有赋值，则实例化对象的该属性是double空矩阵。
disp(class(a));

%%
%对象一旦实例化，其独立属性就生成了，固定不变，除非使用dot运算符或者set方法人为改变。
%从属属性每次被访问都要重新计算，所以如果计算过程中含有随机数，则即使是同一个对象，每次访问结果都很可能不同。
A = [];
B = [];
for i = 1:10
    A = [A;obj.value1];
    B = [B;obj.value2];
end
disp([A,B]);

%%
fprintf('正态分布的标准差为：%d\n',N_randomNumber.sigma);   %Constant属性可以通过类直接调用。
fprintf('正态分布的方差为:%d\n',N_randomNumber.variance);   %一个Constant属性可以从属于另一个Constant属性。

