clc;
clear all;
[vers1, tris1] = readOBJ('./data/spot.obj');
[vers2, tris2] = readOBJ('./data/spot_low_resolution.obj');



%% 设置渲染参数：
shadingParams = {'FaceLighting','gouraud', 'FaceColor','interp'};

% 令每个顶点颜色值正比于其z坐标；
zmax = max(vers1(:, 3));
zmin = min(vers1(:, 3));
range = zmax - zmin;
colorValue = vers1(:, 3)-zmin;
colorValue = colorValue/range;


% tsurf(tris, vers, ....)——三维网格绘图接口
t = tsurf(tris1, vers1, 'CData', colorValue, shadingParams{:});
shading interp;
axis equal;


%% 设置光照

lights = camlight;              % 相机镜头发出的光照

%  light(pram1, value1, param2, value2 ...)
%   Position——光源在三维空间中的位置
%   style
%               local——局部光照
% light('Position',[-1.5 1 1],'Style','local');           % 设置光源属性


%% 颜色设置
colormap(cbrewer('Blues', 500));




