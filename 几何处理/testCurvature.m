clc;
clear all;
close all;

%%
[vers, tris] = readOBJ('./data/��������.obj');
curvs = discrete_mean_curvature(vers, tris);
[adjAngles, ~] = adjacency_dihedral_angle_matrix(vers, tris);


