clc
clear all
close all
dbstop if error

name = 'inputData/model2';
[vers, tris] = readOBJfast([name, '.obj']);
writeOBJ('inputMesh.obj', vers, tris);

%% 调用SDFGen程序生成符号距离场 
tic
command = ['SDFGen.exe', ' ', name, '.obj', ' ', '1.0 3'];
[status, result] = system( command );
fprintf('SDF calculation takes %f s time.\n', toc);


%% 提取SDFGen结果数据
sdfHandle = fopen([name, '.sdf']);              % SDFGen输出数据文件的句柄；
dim = fscanf( sdfHandle, '%d', 3 );         % XYZ三个方向上方块的数量；
ori = fscanf( sdfHandle, '%f', 3 );            % 栅格原点；
space = fscanf( sdfHandle, '%f', 1 );       % 栅格中单个方块的尺寸；
sdf = fscanf( sdfHandle, '%f' );
sdf = reshape(sdf, dim');

% for debug
gridCenters = zeros(dim');
oriIdx = [];

sdf = imfilter(sdf, fspecial('gaussian',3,0.5), 'replicate');
fclose(sdfHandle);


%% 

% 生成三维空间划分网格[X, Y, Z]
[X,Y,Z] = meshgrid(ori(2):space:ori(2)+space*(dim(2)-1), ...
    ori(1):space:ori(1)+space*(dim(1)-1), ...
    ori(3):space:ori(3)+space*(dim(3)-1) );

tic
[trisOut, versOut] = isosurface(X, Y, Z, sdf, -1);
versOut = versOut(:,[2,1,3]);
fprintf('Isosurface extraction takes %f s time.\n', toc);

if(0)
    figure
    drawMesh(V, F, 'facecolor','y', 'edgecolor','none', 'facealpha',1.0);
    drawMesh(Vout, Fout, 'facecolor','g', 'edgecolor','none', 'facealpha',0.5);

    view(3)
    axis equal
    % axis off
    camlight
    lighting gouraud
    cameramenu
    set(gca, 'Position', [0 0 1 1]);
else
	trisOut = trisOut(:, [1, 3, 2]);            % 翻转三角片
    writeOBJ('model_out.obj', versOut, trisOut);
end


%% by Tora——使用一个平面切割网格底部，然后补上需要连接的部分，得到抽壳网格


disp('finished');
