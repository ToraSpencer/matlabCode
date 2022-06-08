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
cubeMesh = struct('vertex',[],'face',[]);

%   cube()——生成立方体网格；
[cubeMesh.vertex, cubeMesh.face] = cube(3, 4, 5);
writeOBJ('cube.obj', cubeMesh.vertex, cubeMesh.face);

 
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

