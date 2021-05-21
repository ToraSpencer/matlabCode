
%���Է����鶼����д��Ax = b��AΪϵ�����󣬶����������B = [A,b]
% �޽� ? R(A)<R(B)
% ��Ψһ�� ? R(A)=R(B)=n
% �������� ? R(A)=R(B)<n
clc;
clear all;
A = [1 2 0 3;3 3 -3 2;2 -2 -3 2;-2 3 2 0];
b = [3 3 5 3]';
B = [A,b];


%%
%���ж��Ƿ��н⡣
%���ǲ�����η��̣�����ǣ�ֱ����A����ռ�
%���������η��̣��ж�����Ψһ�⻹�������⡣��Ψһ��Ļ���ֱ����linsolve���������
%����������ķ���η����飬������ؽ�p���������Ӧ����η�����Ľ�ռ�ı�׼��������
%��η�����Ľ⼯������pƽ�Ƶõ�����η�����Ľ⼯��
%�ؽ�pΪ�������B����Gauss-Jordan��Ԫ�����ɵļ��н��ݾ�������һ�С�
if rank(A)<rank(B)
    disp('�������޽�');
elseif sum(abs(b))==0
        solutionSpace = null(A);
        disp('���������������⣬��ռ�ı�׼������Ϊ��');
        disp(solutionSpace);
elseif rank(A)==size(A,2)
        x = linsolve(A,b);
        disp('��������Ψһ�⣬��Ϊ��');
        disp(x);
else
        solutionSpace = null(A);
        B_prime = rref(B);          %reduced row echelon form
        p = B_prime(:,size(B_prime,2));
        disp('�÷�������Է�������ؽ�Ϊ��');
        disp(p);
        disp('��Ӧ����η�����Ľ�ռ�ı�׼������Ϊ��');
        disp(solutionSpace);
end

        
        
    
        


        

