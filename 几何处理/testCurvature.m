clc;
clear all;
close all;

%%
[vers, tris] = readOBJ('./data/bunny.obj');
curvs = discrete_mean_curvature(vers, tris);
[adjAngles, ~] = adjacency_dihedral_angle_matrix(vers, tris);
shadingParams = {'FaceLighting','gouraud', 'FaceColor','interp'};

% ��ÿ��������ɫֵ�����������ʣ�
cmax = max(curvs);
cmin = min(curvs);
range = cmax - cmin;
colorValue = curvs-cmin;
colorValue = colorValue/range;


% tsurf(tris, vers, ....)������ά�����ͼ�ӿ�
figure(1)
t = tsurf(tris, vers, 'CData', colorValue, shadingParams{:});
shading interp;
axis equal;

light('Position',[-1.5 1 1],'Style','ambient');           % ���ù�Դ����
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
light('Position',[-1.5 1 1],'Style','ambient');           % ���ù�Դ����
colormap(cbrewer('Blues', 500));

