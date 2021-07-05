clc;
clear all;
load('A.mat');
load('B.mat');
load('Acon.mat');
load('Bcon.mat');

temp = sparse(1:length(Acon), Acon, 1, length(Acon), size(A, 2));
tempElems = nonzeros(temp);

A(Acon,:) = temp;
B(Acon,:) = Bcon;

Aelems = nonzeros(A);

b = B(:, 1);
x = A\b;
