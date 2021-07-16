clc;
clear all;
close all;

%%
[V1,F1] = subdivided_sphere(4);
[Vcube,Fcube] = cube(2,2,2);
 Vcube = 1.6.*Vcube;
Vcube = Vcube - 0.8 *ones(size(Vcube));

[W,H] = mesh_boolean(V1,F1,Vcube,Fcube,'intersection');
tsurf(H,W); axis equal;