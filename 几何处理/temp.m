clc;
clear all;
close all;

%%
[V,F,Q] = cube(2, 2, 2);
writeOBJ('cube.obj', V, F);

disp('finished.');