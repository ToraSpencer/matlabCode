clc;
clear all;
%%
%对类进行实例化，产生一个对象。
vector0 = classSample_vector2D();
disp('调用default constructor产生的二维向量对象');
disp(vector0);

vector1 = classSample_vector2D(1,2);
disp('输入参数调用构造方法产生的二维向量对象');
disp(vector1);


vector1.normalize();
disp('调用归一化方法');
disp(vector1);

vector2 = classSample_vector3D(1,2,3);
disp('输入参数调用构造方法产生的三维向量对象：');
disp(vector2);