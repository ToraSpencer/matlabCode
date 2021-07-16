clc;
close all;
clear all;
%%
[vers, tris] = readOBJ('./data/bunny.obj');
oriVers = vers;
oriTris = tris;

L = -cotmatrix(vers,tris);
M = massmatrix(vers,tris);

L = -cotmatrix(vers,tris);
M = massmatrix(vers,tris);
u = 1e-2*L + M;

u = ones(size(vers, 1), 1);

figure(1);
shadingParams = {'FaceLighting','gouraud', 'FaceColor','interp'};
t1 = tsurf(oriTris, oriVers, 'CData',u, shadingParams{:});
shading interp;
axis equal;
lights = camlight;  
colormap(cbrewer('Blues', 500));

figure(2);
t2 = tsurf(tris, vers, 'CData',u, shadingParams{:});
shading interp;
axis equal;
lights = camlight;  
colormap(cbrewer('Blues', 500));

