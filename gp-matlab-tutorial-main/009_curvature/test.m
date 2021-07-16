clc;
clear all;
close all;
%%
[vers, tris] = readOBJ('./data/armadillo.obj');
struct = draw_with_grime(vers, tris);