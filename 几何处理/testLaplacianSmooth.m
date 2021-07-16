clc;
clear all;
close all;

%%
[vers, tris] = readOBJ('./data/bunny.obj');
versNorm = compute_vertex_normals(vers, tris);
L = -cotmatrix(vers,tris);
M = massmatrix(vers,tris);
Melem = nonzeros(M);
Mdiag = diag(M);

temp = L(1, :);
 
b = [];
lambda = 0.1;
[vers_new, vers_new_all] = laplacian_smooth(vers, tris, 'cotan', b,  lambda, 'implicit', vers, 20);
writeOBJ('laplacian_smooth_iter_5.obj', vers_new_all(:, :, 5), tris);
writeOBJ('laplacian_smooth_iter_10.obj', vers_new_all(:, :, 10), tris);
writeOBJ('laplacian_smooth_iter_15.obj', vers_new_all(:, :, 15), tris);
writeOBJ('laplacian_smooth_iter_20.obj', vers_new, tris);

