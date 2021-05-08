clc
clear all;
A = rand(10,1);
B = A;
C = 1:10;
C = C';
tableHead = {'a','b','c'};
T = table(A,B,C,'VariableNames',tableHead);
disp(T);

%%
T_inverted = sortrows(T,{'c'},{'descend'});
A = T_inverted.a;
B = T_inverted.b;
C = T_inverted.c;
disp(T_inverted);