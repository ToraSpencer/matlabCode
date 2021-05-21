
%线性方程组都可以写作Ax = b，A为系数矩阵，定义增广矩阵B = [A,b]
% 无解 ? R(A)<R(B)
% 有唯一解 ? R(A)=R(B)=n
% 有无穷多解 ? R(A)=R(B)<n
clc;
clear all;
A = [1 2 0 3;3 3 -3 2;2 -2 -3 2;-2 3 2 0];
b = [3 3 5 3]';
B = [A,b];


%%
%先判断是否有解。
%再是不是齐次方程，如果是，直接求A的零空间
%如果不是齐次方程，判断是有唯一解还是无穷多解。有唯一解的话，直接用linsolve函数解出。
%有无穷多个解的非齐次方程组，先求出特解p，再求出对应的齐次方程组的解空间的标准正交基。
%齐次方程组的解集经向量p平移得到非齐次方程组的解集。
%特解p为增广矩阵B经过Gauss-Jordan消元法化成的简化行阶梯矩阵的最后一列。
if rank(A)<rank(B)
    disp('方程组无解');
elseif sum(abs(b))==0
        solutionSpace = null(A);
        disp('方程组有无穷多个解，解空间的标准正交基为：');
        disp(solutionSpace);
elseif rank(A)==size(A,2)
        x = linsolve(A,b);
        disp('方程组有唯一解，解为：');
        disp(x);
else
        solutionSpace = null(A);
        B_prime = rref(B);          %reduced row echelon form
        p = B_prime(:,size(B_prime,2));
        disp('该非齐次线性方程组的特解为：');
        disp(p);
        disp('对应的齐次方程组的解空间的标准正交基为：');
        disp(solutionSpace);
end

        
        
    
        


        

