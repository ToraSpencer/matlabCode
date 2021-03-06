clc
clear all
close all
dbstop if error


%%
type = 1;           % 1:Union 2:Intersection 3:Difference
path = 'myData\';  
file = {'', '', ''};
 
file(1) = {[path, 's9block_preprocessed.obj']};
file(2) = {[path, 'fatTeeth1_preprocessed.obj']};
file(3) = {[path, 'fatTeeth2_preprocessed.obj']};
meshes = cell(2,3);

for i = 1:3
    [meshes{1,i}, meshes{2,i}] = readOBJ(file{i});
end


%%
% 标记
nfile = 3;
V0 = meshes{1,1};
F0 = meshes{2,1};
Tag = zeros(size(F0,1),1);

for j = 2:nfile
    V0 = V0(F0(:),:);
    F0 = reshape(1:3*size(F0,1), size(F0,1), 3);

    V1 = meshes{1,j};
    F1 = meshes{2,j};
    
    V1 = V1(F1(:),:);
    F1 = reshape(1:3*size(F1,1), size(F1,1), 3);

    Tag = [Tag; (j-1)*ones(size(F1,1),1)];
    F0 = [F0; F1 + size(V0,1)];
    V0 = [V0; V1];
end


%% 第三方库程序meshArrange处理——三角剖分：
writeSTL('input.stl', V0, F0);
dlmwrite('inputlabel.txt', Tag, 'delimiter','\t', 'newline','pc', 'precision','%.0f');
tic
[status,result] = system('runme.bat');
fprintf('mesh arrangement takes %f s time.\n', toc);

[vers, tris] = readOBJfast('output.obj');
outputlabel = load('outputlabel.txt');
outputlabel = log10(outputlabel) + 1;


%% 处理自相交——舍弃内部三角片，只去外层三角片，获得union结果；
if type == 1
    tic
    
    temp = tris(:, 1);
    temp = vers(temp,3);
    [~, fidx] = min(temp);          % 网格中z坐标最小的顶点所在的三角片索引，貌似是刻意选一个交叉部分以外的三角片；
    theTri = tris(fidx, :);
    objWriteVertices('chosenVer.obj', vers(theTri', :));
    [Vout, Fout] = choose_tris(vers, tris, fidx);
    
%     [Vout, Fout] = meshfix(Vout, Fout);
    fprintf('solve self intersection takes %f s time.\n', toc);

    writeOBJ('finalMesh.obj', Vout, Fout);
end

figure
drawMesh(Vout, Fout, 'facecolor','y', 'edgecolor','none', 'facealpha',0.9);

view(3)
axis equal
axis off
camlight
lighting gouraud
cameramenu
set(gca, 'Position', [0 0 1 1]);
