clc
clear all
close all
dbstop if error

name = 'inputData/model2';
[vers, tris] = readOBJfast([name, '.obj']);
writeOBJ('inputMesh.obj', vers, tris);

%% ����SDFGen�������ɷ��ž��볡 
tic
command = ['SDFGen.exe', ' ', name, '.obj', ' ', '1.0 3'];
[status, result] = system( command );
fprintf('SDF calculation takes %f s time.\n', toc);


%% ��ȡSDFGen�������
sdfHandle = fopen([name, '.sdf']);              % SDFGen��������ļ��ľ����
dim = fscanf( sdfHandle, '%d', 3 );         % XYZ���������Ϸ����������
ori = fscanf( sdfHandle, '%f', 3 );            % դ��ԭ�㣻
space = fscanf( sdfHandle, '%f', 1 );       % դ���е�������ĳߴ磻
sdf = fscanf( sdfHandle, '%f' );
sdf = reshape(sdf, dim');

% for debug
gridCenters = zeros(dim');
oriIdx = [];

sdf = imfilter(sdf, fspecial('gaussian',3,0.5), 'replicate');
fclose(sdfHandle);


%% 

% ������ά�ռ仮������[X, Y, Z]
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
	trisOut = trisOut(:, [1, 3, 2]);            % ��ת����Ƭ
    writeOBJ('model_out.obj', versOut, trisOut);
end


%% by Tora����ʹ��һ��ƽ���и�����ײ���Ȼ������Ҫ���ӵĲ��֣��õ��������


disp('finished');
