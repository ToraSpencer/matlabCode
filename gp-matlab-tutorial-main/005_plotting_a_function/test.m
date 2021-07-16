clc;
clear all;
[vers, tris] = readOBJ('./data/spot.obj');
plot_z_coord(vers, tris);