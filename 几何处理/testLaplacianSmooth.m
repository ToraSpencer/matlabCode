clc;
clear all;
close all;
[vers, tris] = readOBJ('./data/bunny.obj');

% ∆Ωª¨

%% laplacian
clc;
clear all;
close all;
[vers, tris] = readOBJ('./data/bunny.obj');
L = cotmatrix_embedded(vers,tris);


%% ≤‚ ‘laplace∆Ωª¨
clc;
clear all;
close all;
[vers, tris] = readOBJ('./data/bunny.obj');
versNorm = compute_vertex_normals(vers, tris);

b = [];
lambda = 0.1;

% laplacian_smooth()
[vers_new, vers_new_all] = laplacian_smooth(vers, tris, 'cotan', b,  lambda, 'implicit', vers, 20);
writeOBJ('laplacian_smooth_iter_5.obj', vers_new_all(:, :, 5), tris);
writeOBJ('laplacian_smooth_iter_10.obj', vers_new_all(:, :, 10), tris);
writeOBJ('laplacian_smooth_iter_15.obj', vers_new_all(:, :, 15), tris);
writeOBJ('laplacian_smooth_iter_20.obj', vers_new, tris);




