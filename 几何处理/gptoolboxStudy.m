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
cubeMesh = struct('vertex',[],'face',[]);

%   cube()������������������
[cubeMesh.vertex, cubeMesh.face] = cube(3, 4, 5);
writeOBJ('cube.obj', cubeMesh.vertex, cubeMesh.face);

 
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

