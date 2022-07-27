%% gptoolbox - Geometry Processing Toolbox
 clc;
 clear all;

 
%% ��ܵ�����ƽ���㷨��
clc;
clear all;
beforeSmooth = readOBJ('beforeSmooth.obj');
param = 0.01;
%        smooth_loop()�������޶���ά�ռ䣬����ά��Ҳ���ԡ�
versCount = size(beforeSmooth,1);
vec1 = mod((1:versCount)-2, versCount) + 1;                 % mod()ȡ����
vec2 = mod(1:versCount, versCount) + 1;
vec3 = [vec1, vec2];
W = sparse([1:versCount, 1:versCount], vec3, 1);

A1 = speye(versCount) - 0.5 * W;
B1 = zeros(size(beforeSmooth));
A2 = speye(versCount);
B2 = beforeSmooth;
A = [A1; param*A2];         %2v*v
B = [B1; param*B2];         %2v*3
afterSmooth = A\B;

temp1 = A1\B1;
temp2 = A2\B2;
 
OBJwriteVertices('afterSmooth2.obj', afterSmooth);



%% IO�ӿڣ�
tooth =  Read_Obj('./data/tooth.obj');      

% writeOBJ()����дOBJ�ļ���
writeOBJ('newTooth.obj', tooth.vertex, tooth.face);
objWriteVertices('toothVers.obj', tooth.vertex);
 


%% ���ɻ���ͼ�Σ�
clc;
clear all;
cubeMesh = struct('vertex',[],'face',[]);

%   cube()������������������
[cubeMesh.vertex, cubeMesh.face] = cube(3, 4, 5);           % ��(0,0,0)��ʼ����xyz���궼���ڵ���0
writeOBJ('cube.obj', cubeMesh.vertex, cubeMesh.face);

%    interpolateToLine()������ֵ����ֱ�߶ε㼯����by Tora
line = interpolateToline([0, 0, 0], [3, 4, 5], 0.5);
objWriteVertices('line.obj', line);

%    printAABB()�������ɲ�ΰ�Χ������д�뵽�����ļ�����by Tora
printAABB('aabb.obj', [-4.65682888, -6.20851517, -5.11108398], [4.63390493, 4.79591894, 3.49578595]);

 
%% ��ͼ
clc;
clear all;
tooth =  Read_Obj('./data/tooth.obj');  
toothVers = tooth.vertex;

%       ����ά����
figure(1);
trimesh(tooth.face,tooth.vertex(:,1),tooth.vertex(:,2),tooth.vertex(:,3));

%       ������
figure(2)
plot3(toothVers(:, 1), toothVers(:, 2), toothVers(:, 3), '.');

%        getDirLine()������ȡ����ָʾ�ߵ㼯
xline = getDirLine([1, 0, 0], [0, 0, 0]);
figure(3);
plot3(xline(:, 1), xline(:, 2), xline(:, 3), '*');



%% MyCrustOpen()����������ƣ��������Ƭ����������ò�����ɵ�����Ƭ�������
trimesh(MyCrustOpen(toothVers), toothVers(:, 1), toothVers(:, 2), toothVers(:, 3));
writeOBJ('������Ƭ���ɵ�tooth����.obj', toothVers, MyCrustOpen(toothVers));



%% �鶴������
clc;
clear all;
[vers, tris] = readOBJ('�и�����.obj');
bdry =  detect_mesh_holes_and_boundary(tris);
holeVersIdx = bdry{1,1};
holeVers = vers(holeVersIdx, :);
objWriteVertices('�и����ڵĶ�.obj', holeVers);



%% ��������Ļ������ʣ�
clc;
clear all;
[vers, tris] = readOBJfast('data/spot_low_resolution.obj');
writeOBJ('inputMesh.obj', vers, tris);
versCount = size(vers, 1);
trisCount = size(tris, 1);

%                            ȷ�����ʵ�ָʾ�߳��ȣ�
minX = min(vers(:, 1));
maxX = max(vers(:, 1));
minY = min(vers(:, 2));
maxY = max(vers(:, 2));
minZ = min(vers(:, 3));
maxZ = max(vers(:, 3));
lineLen = max([maxX - minX, maxY - minY, maxZ - minZ]);

% barycenter()����������������ÿ������Ƭ�����ģ�
baryCenters = barycenter(vers, tris);

% normals()��������������Ƭ�ķ���
normDirs = normals(vers, tris, 'UseSVD', 1);                % δ��һ����
normDirs = normalizerow(normDirs);
%                           �����㷨��ָʾ�ߣ�
startVers = baryCenters;
endVers = startVers + normDirs * lineLen;
lines = [];
for i = 1: trisCount
    currentLine = interpolateToline(startVers(i, :), endVers(i, :), lineLen/10);
    lines = [lines; currentLine];
end
objWriteVertices('trisNormDirLines.obj', lines);
objWriteVertices('firstTrisNormDirLine.obj', lines(1: floor(size(lines, 1)/trisCount), :));

% per_vertex_normals()���������񶥵㷨��
normDirs_ver = per_vertex_normals(vers, tris);
%                       �����㷨��ָʾ�ߣ�
startVers = vers;
endVers = startVers + normDirs_ver*lineLen;
lines = [];
for i = 1: versCount
    currentLine = interpolateToline(startVers(i, :), endVers(i, :), lineLen/10);
    lines = [lines; currentLine];
end
objWriteVertices('versNormDirLines.obj', lines);
objWriteVertices('firstVersNormDirLine.obj', lines(1: floor(size(lines, 1)/versCount), :));




%% ��������Ļ������ʣ�
clc;
clear all;


%%
 clc;
 clear all;

[cubeMesh.vertex, cubeMesh.face] = cube(3, 4, 5);
writeOBJ('cube.obj', cubeMesh.vertex, cubeMesh.face);

edges = [1,2; 2,3; 3,4; 8,9; 9,10; 11,12];
objWriteEdges('cubeEdges.obj', edges, cubeMesh.vertex);

ret = connected_components(cubeMesh.face);

