%% gptoolbox - Geometry Processing Toolbox
 

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
 
