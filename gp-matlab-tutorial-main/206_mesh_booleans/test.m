clc;
clear all;
close all;

[vers1,tris1] = subdivided_sphere(4);
writeOBJ('./sphere.obj', vers1, tris1);

[Vcube,Fcube] = cube(2,2,2);
Vcube = 1.6.*Vcube;
Vcube = Vcube - 0.8 *ones(size(Vcube));
writeOBJ('./cube.obj', Vcube, Fcube);

shadingParams = {'FaceLighting','gouraud', 'FaceColor','interp'};

% ��ÿ��������ɫֵ��������z���ꣻ
zmax = max(vers1(:, 3));
zmin = min(vers1(:, 3));
range = zmax - zmin;
[rows, cols] = size(vers1);
colorValue = ones(rows, 1);

colorValue = colorValue/range;
t = tsurf(tris1, vers1, 'CData', colorValue, shadingParams{:});
shading interp;
axis equal;
lights = camlight;              % �����ͷ�����Ĺ���


%%
[W,H] = mesh_boolean(vers1,tris1,Vcube,Fcube,'union');
tsurf(H,W); axis equal;