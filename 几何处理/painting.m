clc;
clear all;
[vers1, tris1] = readOBJ('./data/bunny.obj');


%% ��������matlab�Դ���plot3()�ӿڣ�
plot3(vers1(:, 1), vers1(:, 2), vers1(:, 3), '.');
 


%% ������Ⱦ������
shadingParams = {'FaceLighting','gouraud', 'FaceColor','interp'};

% ��ÿ��������ɫֵ��������z���ꣻ
zmax = max(vers1(:, 3));
zmin = min(vers1(:, 3));
range = zmax - zmin;
colorValue = vers1(:, 3)-zmin;
colorValue = colorValue/range;


% tsurf(tris, vers, ....)������ά�����ͼ�ӿ�
t = tsurf(tris1, vers1, 'CData', colorValue, shadingParams{:});
shading interp;
axis equal;


%% ���ù���

lights = camlight;              % �����ͷ�����Ĺ���

%  light(pram1, value1, param2, value2 ...)
%   Position������Դ����ά�ռ��е�λ��
%   style
%               local�����ֲ�����
% light('Position',[-1.5 1 1],'Style','local');           % ���ù�Դ����


%% ��ɫ����
% colormap(cbrewer('Blues', 500));




