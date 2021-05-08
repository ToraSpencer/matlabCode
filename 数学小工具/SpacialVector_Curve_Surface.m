%%
%��ά�ռ�����������
clc;
clear all;
syms x y z t
u = x^2+log(y)+z^2 ;
v = [x^2+y+z y+log(x) z+2*x*y];
gradient = [diff(u,'x') diff(u,'y') diff(u,'z')];
rotation = [diff(v(3),'y')-diff(v(2),'z'),diff(v(1),'z')-diff(v(3),'x'),diff(v(2),'x')-diff(v(1),'y')];
divergence = diff(v(1),'x')+diff(v(2),'y')+diff(v(3),'z');
disp('��������uΪ��');
disp(u);
disp('u���ݶ�Ϊ:');
disp(gradient);
disp('��������vΪ��');
disp(v);
disp('v������Ϊ��');
disp(rotation);
disp('v��ɢ��Ϊ��');
disp(divergence);


%%
%��ά�ռ�����������:
clc;
clear all;
a = [1,2,3];
b = [3,4,5];
innerProduct = dot(a,b);
outerProduct = cross(a,b);
disp(sprintf('�����������ڻ�Ϊ��%g\n�������������Ϊ��(%g,%g,%g)',innerProduct,outerProduct));

%%
%�ռ�������ز������㣺
clc;
clear all;
syms t;
r = [3*sin(t) 3*cos(t) -2*t];
v = diff(r);
T = v/norm(v);              %���ߵĵ�λ��������
disp('T = ');
disp(T);

kappa = 1/norm(v)*norm(diff(T));        %���ߵ����ʡ�
disp('kappa = ');
disp(kappa);

rou = 1/kappa;              %���ߵ����ʰ뾶��
disp('rou = ');
disp(rou);

N = diff(T)/norm(diff(T));  %���ߵ�����λ��������
disp('N = ');
disp(N);

B = cross(T,N);             %���ߵĸ���������
disp('B = ');
disp(B);
