clc
clear all
close all
dbstop if error


%%
type = 1;           % 1:Union 2:Intersection 3:Difference
path = 'data\data1\'; % data1: 原型
file = dir([path, '*.obj']);
nfile = length(file);
model = cell(2,nfile);

temp = [path, file(1).name];
for i = 1:nfile
    [model{1,i}, model{2,i}] = readOBJ([path, file(i).name]);
end


%%
% 布尔并操作；
V0 = model{1,1};
F0 = model{2,1};
Tag = zeros(size(F0,1),1);

for j = 2:nfile
    V0 = V0(F0(:),:);
    F0 = reshape(1:3*size(F0,1), size(F0,1), 3);

    V1 = model{1,j};
    F1 = model{2,j};
    
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


%% 最终结果输出
% 读取mesh arrangement结果
[V, F] = readOBJfast('output.obj');
outputlabel = load('outputlabel.txt');
outputlabel = log10(outputlabel) + 1;


if type == 1
    tic
    [useless, fidx] = min(V(F(:,1),3));
    [Vout, Fout] = choose_tris(V, F, fidx);
    
%     [Vout, Fout] = meshfix(Vout, Fout);
    fprintf('solve self intersection takes %f s time.\n', toc);

    writeOBJ([path, 'result\model_out.obj'], Vout, Fout);
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
