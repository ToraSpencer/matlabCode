clc;
clear all;
close all;

%%
[vers, tris] = readOBJ('./data/bunny.obj');
curvs = discrete_mean_curvature(vers, tris);
[adjAngles, ~] = adjacency_dihedral_angle_matrix(vers, tris);
shadingParams = {'FaceLighting','gouraud', 'FaceColor','interp'};

% 令每个顶点颜色值正比于其曲率；
cmax = max(curvs);
cmin = min(curvs);
range = cmax - cmin;
colorValue = curvs-cmin;
colorValue = colorValue/range;


% tsurf(tris, vers, ....)——三维网格绘图接口
figure(1)
t = tsurf(tris, vers, 'CData', colorValue, shadingParams{:});
shading interp;
axis equal;

light('Position',[-1.5 1 1],'Style','ambient');           % 设置光源属性
colormap(cbrewer('Blues', 500));


b = [];
lambda = 0.1;
[vers_new, vers_new_all] = laplacian_smooth(vers, tris, 'cotan', b,  lambda, 'implicit', vers, 5);
tris_new = tris;
curvs = discrete_mean_curvature(vers_new, tris_new);
shadingParams = {'FaceLighting','gouraud', 'FaceColor','interp'};
cmax = max(curvs);
cmin = min(curvs);
range = cmax - cmin;
colorValue = curvs - cmin;
colorValue = colorValue/range;
figure(2)
t = tsurf(tris_new, vers_new, 'CData', colorValue, shadingParams{:});
shading interp;
axis equal;
light('Position',[-1.5 1 1],'Style','ambient');           % 设置光源属性
colormap(cbrewer('Blues', 500));

