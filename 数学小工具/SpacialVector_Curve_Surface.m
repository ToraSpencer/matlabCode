%%
%三维空间向量分析：
clc;
clear all;
syms x y z t
u = x^2+log(y)+z^2 ;
v = [x^2+y+z y+log(x) z+2*x*y];
gradient = [diff(u,'x') diff(u,'y') diff(u,'z')];
rotation = [diff(v(3),'y')-diff(v(2),'z'),diff(v(1),'z')-diff(v(3),'x'),diff(v(2),'x')-diff(v(1),'y')];
divergence = diff(v(1),'x')+diff(v(2),'y')+diff(v(3),'z');
disp('标量函数u为：');
disp(u);
disp('u的梯度为:');
disp(gradient);
disp('向量函数v为：');
disp(v);
disp('v的旋度为：');
disp(rotation);
disp('v的散度为：');
disp(divergence);


%%
%三维空间向量的运算:
clc;
clear all;
a = [1,2,3];
b = [3,4,5];
innerProduct = dot(a,b);
outerProduct = cross(a,b);
disp(sprintf('两个向量的内积为：%g\n两个向量的外积为：(%g,%g,%g)',innerProduct,outerProduct));

%%
%空间曲线相关参数运算：
clc;
clear all;
syms t;
r = [3*sin(t) 3*cos(t) -2*t];
v = diff(r);
T = v/norm(v);              %曲线的单位切向量。
disp('T = ');
disp(T);

kappa = 1/norm(v)*norm(diff(T));        %曲线的曲率。
disp('kappa = ');
disp(kappa);

rou = 1/kappa;              %曲线的曲率半径。
disp('rou = ');
disp(rou);

N = diff(T)/norm(diff(T));  %曲线的主单位法向量。
disp('N = ');
disp(N);

B = cross(T,N);             %曲线的副法向量。
disp('B = ');
disp(B);
