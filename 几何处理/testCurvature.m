clc;
clear all;
close all;

%%
[vers, tris] = readOBJ('./data/正四面体.obj');
curvs = discrete_mean_curvature(vers, tris);
[adjAngles, ~] = adjacency_dihedral_angle_matrix(vers, tris);


