%% gptoolbox - Geometry Processing Toolbox
 clc;
 clear all;

 
%% 斌杰的曲线平滑算法：
clc;
clear all;
beforeSmooth = readOBJ('beforeSmooth.obj');
param = 0.01;
%        smooth_loop()——不限定三维空间，其他维数也可以。
versCount = size(beforeSmooth,1);
vec1 = mod((1:versCount)-2, versCount) + 1;                 % mod()取余数
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



%% IO接口：
tooth =  Read_Obj('./data/tooth.obj');      

% writeOBJ()——写OBJ文件；
writeOBJ('newTooth.obj', tooth.vertex, tooth.face);
objWriteVertices('toothVers.obj', tooth.vertex);
 


%% 生成基础图形：
clc;
clear all;
cubeMesh = struct('vertex',[],'face',[]);

%   cube()——生成立方体网格；
[cubeMesh.vertex, cubeMesh.face] = cube(3, 4, 5);           % 从(0,0,0)开始画，xyz坐标都大于等于0
writeOBJ('cube.obj', cubeMesh.vertex, cubeMesh.face);

%    interpolateToLine()——插值生成直线段点集——by Tora
line = interpolateToline([0, 0, 0], [3, 4, 5], 0.5);
objWriteVertices('line.obj', line);

%    printAABB()——生成层次包围盒网格，写入到本地文件——by Tora
printAABB('aabb.obj', [-4.65682888, -6.20851517, -5.11108398], [4.63390493, 4.79591894, 3.49578595]);

 
%% 画图
clc;
clear all;
tooth =  Read_Obj('./data/tooth.obj');  
toothVers = tooth.vertex;

%       画三维网格
figure(1);
trimesh(tooth.face,tooth.vertex(:,1),tooth.vertex(:,2),tooth.vertex(:,3));

%       画点云
figure(2)
plot3(toothVers(:, 1), toothVers(:, 2), toothVers(:, 3), '.');

%        getDirLine()——获取方向指示线点集
xline = getDirLine([1, 0, 0], [0, 0, 0]);
figure(3);
plot3(xline(:, 1), xline(:, 2), xline(:, 3), '*');



%% MyCrustOpen()——输出点云，输出三角片——！！！貌似生成的三角片朝向混乱
trimesh(MyCrustOpen(toothVers), toothVers(:, 1), toothVers(:, 2), toothVers(:, 3));
writeOBJ('补三角片生成的tooth网格.obj', toothVers, MyCrustOpen(toothVers));



%% 查洞、补洞
clc;
clear all;
[vers, tris] = readOBJ('切割牙冠.obj');
bdry =  detect_mesh_holes_and_boundary(tris);
holeVersIdx = bdry{1,1};
holeVers = vers(holeVersIdx, :);
objWriteVertices('切割牙冠的洞.obj', holeVers);



%% 三角网格的基本性质：
clc;
clear all;
[vers, tris] = readOBJfast('data/spot_low_resolution.obj');
writeOBJ('inputMesh.obj', vers, tris);
versCount = size(vers, 1);
trisCount = size(tris, 1);

%                            确定合适的指示线长度：
minX = min(vers(:, 1));
maxX = max(vers(:, 1));
minY = min(vers(:, 2));
maxY = max(vers(:, 2));
minZ = min(vers(:, 3));
maxZ = max(vers(:, 3));
lineLen = max([maxX - minX, maxY - minY, maxZ - minZ]);

% barycenter()——计算三角网格每个三角片的重心：
baryCenters = barycenter(vers, tris);

% normals()——计算网格面片的法向
normDirs = normals(vers, tris, 'UseSVD', 1);                % 未归一化；
normDirs = normalizerow(normDirs);
%                           画顶点法向指示线：
startVers = baryCenters;
endVers = startVers + normDirs * lineLen;
lines = [];
for i = 1: trisCount
    currentLine = interpolateToline(startVers(i, :), endVers(i, :), lineLen/10);
    lines = [lines; currentLine];
end
objWriteVertices('trisNormDirLines.obj', lines);
objWriteVertices('firstTrisNormDirLine.obj', lines(1: floor(size(lines, 1)/trisCount), :));

% per_vertex_normals()——求网格顶点法向：
normDirs_ver = per_vertex_normals(vers, tris);
%                       画顶点法向指示线：
startVers = vers;
endVers = startVers + normDirs_ver*lineLen;
lines = [];
for i = 1: versCount
    currentLine = interpolateToline(startVers(i, :), endVers(i, :), lineLen/10);
    lines = [lines; currentLine];
end
objWriteVertices('versNormDirLines.obj', lines);
objWriteVertices('firstVersNormDirLine.obj', lines(1: floor(size(lines, 1)/versCount), :));




%% 三角网格的基本性质：
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

