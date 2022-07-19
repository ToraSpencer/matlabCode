clc;
clear all;
close all;
debugFlag = 1;

%% ���������ݣ�prepare data:
load('inputData/data.mat');

aveBone.vers = [];
aveBone.tris = [];

[aveBone.vers, aveBone.tris] = readOBJ('inputData/aveBone_lower_stan.obj');       % ��׼������۹�����
conVers = [                 % ��׼���۹ǵĿ��Ƶ㡣
    -32.0, -36.0, 0
    -28.5, -27.5, 0
    -24.5, -20.0, 0
    -21.0, -13.0, 0
    -15.5, -7.0, 0
    -10.0, -2.0, 0
    -4.0, -0.5, 0
    4.0, -0.5, 0
    10.0, -2.0, 0
    15.5, -7.0, 0
    21.0, -13.0, 0
    24.5, -20.0, 0
    28.5, -27.5, 0
    32.0, -36.0, 0];
conVers = conVers(end:-1:1,:);

dentalCenter = handle.dentalCenter;
dentalFrame = handle.dentalFrame;
toothdata = handle.model.toothdata;

teethCount = size(toothdata,2);

% ��������
teethCenter = zeros(teethCount,3);
for i = 1:teethCount
    teethCenter(i,:) = mean(toothdata{1,i}(toothdata{3,i},:));          % �ֶ����������ߵ������������е�������
end


%% check input
if(debugFlag == 1)
    objWriteVertices('conVers.obj', conVers);
    tooth.vers = toothdata{1, 1};               % ��������
    tooth.tris = toothdata{2, 1};
    tooth.rVers = toothdata{9, 1};              % ������������
    tooth.rTris = toothdata{10, 1};
    writeOBJ('tooth1.obj', tooth.vers, tooth.tris);
    writeOBJ('rootTooth1.obj', tooth.rVers, tooth.rTris);
    objWriteVertices('dentalFrame.obj', dentalFrame);
    objWriteVertices('dentalCenter.obj', dentalCenter);
    objWriteVertices('toothCenter.obj', teethCenter);
end
    

%% ����������ת�������xyƽ���ϲ��ñ�׼��Բ���

%1. �����任�����������ĵ�任���������ϵ��
teethCenterJC = bsxfun(@minus, teethCenter, dentalCenter) * dentalFrame;    

%2. �������ĵ������Բ������׼��Բ��û�������������ת��a*x^2 + c*y^2 + d*x + e*y + f = 0;
sampleVers = teethCenterJC(:, 1:2);
coff = fit_ellipse(sampleVers(:,1), sampleVers(:,2), 'standard');
a = coff(1);
c = coff(2);
d = coff(3);
e = coff(4);
f = coff(5);

center = [-d/(2*a), -e/(2*c)];            % ��Բ�������ꣻ
p = d^2/(4*a) + e^2/(4*c) - f;
a0 = sqrt(p/a);
b0 = sqrt(p/c);

step = (-pi/2:0.01:3*pi/2)';
elliJC = [a0*cos(step), b0*sin(step)];
elliJC = bsxfun(@plus, elliJC, center);

if(0)
    figure
    hold on
    scatter(sampleVers(:,1), sampleVers(:,2), 'r');
    scatter(elliJC(:, 1), elliJC(:, 2), 'b');
end


% 3. �����Բ����������ѡȡ��������֮��ĵ�..
startVer = (sampleVers(1,:) + [-sampleVers(teethCount,1), sampleVers(teethCount,2)])/2;
idx1 = knnsearch(elliJC, startVer);                 % �����Բ�Ͼ���startVer����ĵ��������
idx2 = knnsearch(step, pi-step(idx1));
sel = round(linspace(idx1, idx2, 14));
targetVersJC = elliJC(sel,:);        
objWriteVertices('target1.obj', [elliJC(idx1, :), 0]);
objWriteVertices('target2.obj', [elliJC(idx2, :), 0]);

if(debugFlag == 1)
    objWriteVertices('sampleVers.obj', teethCenterJC);
    objWriteVertices('�����Բ.obj', [elliJC, zeros(size(elliJC, 1), 1)]);
    objWriteVertices('����Ŀ���.obj', [targetVersJC, zeros(size(targetVersJC, 1), 1)]);
end


%% ��׼�Ǳ��Σ�ת����ȫ������ϵ
% ����Ȩ��
temp = reshape(conVers(:,1:2)', [1 2 14]);
weight = bsxfun(@minus, aveBone.vers(:,1:2), temp);
weight = sum(weight.^2, 2)/1000;
weight = exp(-weight);
weight = bsxfun(@rdivide, weight, sum(weight,3));
arrows = (targetVersJC - conVers(:,1:2))';         % ���Ƶ�ָ����Բ������Ķ�ά���������洢�ɾ���
tempMat = reshape(arrows, 1, 2, teethCount);
transMat = bsxfun(@times, tempMat, weight);

% ���Ρ���ֻ��д��׼���۹Ƕ����xy����ֵ��
aveBone.vers(:,1:2) = sum( transMat, 3) + aveBone.vers(:,1:2);

% �任��global����ϵ��
aveBone.vers = bsxfun(@plus, aveBone.vers * dentalFrame', dentalCenter);


%% д������ݣ�
rootTeeth.vers = [];
rootTeeth.tris = [];
 
for i = 1:teethCount
    rootTeeth.tris = [rootTeeth.tris; toothdata{10,i} + size(rootTeeth.vers,1)];
    rootTeeth.vers = [rootTeeth.vers; toothdata{9,i}];
end


if(debugFlag == 1)
    writeOBJ('rootTeethLower.obj', rootTeeth.vers, rootTeeth.tris);
    writeOBJ('finalMesh.obj', aveBone.vers, aveBone.tris);
else
    figure
    drawMesh(rootTeeth.vers, rootTeeth.tris, 'facecolor','y', 'edgecolor','none', 'facealpha', 1.0);
    drawMesh(aveBone.vers, aveBone.tris, 'facecolor','g', 'edgecolor','none', 'facealpha', 0.8);

    view(3)
    axis equal
    axis off
    camlight
    lighting gouraud
    set(gca, 'Position',[0 0 1 1]);
end

disp('main() finished.');
 
 