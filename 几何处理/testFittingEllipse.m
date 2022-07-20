clc;
clear all;
close all;

%%
load('./data/elliSampleVers.mat');
x0 = sampleVers(:, 1);
y0 = sampleVers(:, 2);
sampleCount = size(sampleVers, 1);
figure
scatter(x0, y0, 'r', 'filled');
hold on

%% 最小二乘法拟合椭圆：
%  标准椭圆方程：a*x^2 + c*y^2 + d*x + e*y + f  = 0，其中a*c > 0
%   alpha = [x^2, y^2, x, y, 1]; 样本信息矩阵：A = [alpha1; alpha2; .... alpham]

A = [x0.*x0, y0.*y0, x0, y0, ones(sampleCount, 1)];       
% 标准椭圆方程可写为A*x = 0; x = [a,c,d,e,f];
% 最小二乘法拟合椭圆本质上是在a*c>0的约束下求A*x = 0最优解；
% 法线方程：A'*A*x = 0;
% 余项：r = -A*x;
 
 %%
 B = zeros(5, 5);
 B(1, 2) = 1;
 B(2, 1) = 1;
 S = inv(A'*A)*transpose(B);
 %  S*x = 1/lambda* x; x = [a, c, d, e, f];
 [V, D] = eig(S);
 elemInfo.row = [];
elemInfo.col = [];
elemInfo.value = [];
[elemInfo.row, elemInfo.col, elemInfo.value] = find(D);

% 不为0的特征值可能不止一个，检测是否符合约束条件：
eigenVec = [];
for i = 1 : size(elemInfo.col)
    index = elemInfo.col(i);
    eigenVec = V(:, index);
    a = eigenVec(1);
    c = eigenVec(2);
    if(a*c > 0)
        break;
    end
end

 a = eigenVec(1);
 c = eigenVec(2);
 d = eigenVec(3);
 e = eigenVec(4);
 f = eigenVec(5);
 
 % 椭圆圆心，半长轴半短轴：
 x0 = -d/(2*a);
 y0 = -e/(2*c);
 p = d^2/(4*a) + e^2/(4*c) - f;
 a0 = sqrt(p/a);
 b0 = sqrt(p/c);
 
 % 椭圆的参数方程：x = a0*cos(theta) + x0; y = b0*sin(theta) + y0;
 theta = [0: 0.05: 2*pi];
 x = a0*cos(theta) + x0;
 y = b0*sin(theta) + y0;
 scatter(x, y, 'b');
 
 
 
 
 
 
