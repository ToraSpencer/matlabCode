%%
clc;
close all;
clear all;

x = linspace(-5, 15, 201);

L = 10;

H = 8;
X = x - L/2;
y = 16/(L^4)*H * X.^4 - H;

plot(x, y, 'r');

% plot(x,y);