clc;
clear all;
close all;
debugFlag = 1;

%% ���������ݣ�prepare data:
aveBone.vers = [];
aveBone.tris = [];
conVers = [];
isUpper = 0;

if(1 == isUpper)
     load('inputData/dataUpper.mat');
    [aveBone.vers, aveBone.tris] = readOBJ('inputData/aveBone_upper_stan.obj');       % ��׼������۹�����
    conVers = [                 % ��׼���۹ǵĿ��Ƶ㡣
            -24.0, -30.0, 0
            -22.5, -22.0, 0
            -21.0, -14.5, 0
            -19.0, -8.5, 0
            -15.5, -3.5, 0
            -10.5, -0.5, 0
            -3.5, 1.0, 0
            3.5, 1.0, 0
            10.5, -0.5, 0
            15.5, -3.5, 0
            19.0, -8.5, 0
            21.0, -14.5, 0
            22.5, -22.0, 0
            24.0, -30.0, 0];
    objWriteVertices('conVers.obj', conVers);
else
    load('inputData/data.mat');
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
    objWriteVertices('conVers.obj', conVers);
    conVers = conVers(end:-1:1, :);     % �������������򽫿��Ƶ㵹װ��
end

dentalCenter = handle.dentalCenter;
rotation = handle.dentalFrame;
toothdata = handle.model.toothdata;

teethCount = size(toothdata,2);

% ����������
glCenters = zeros(teethCount,3);
for i = 1:teethCount
    glCenters(i,:) = mean(toothdata{1,i}(toothdata{3,i},:));          % �ֶ����������ߵ������������е�������
end


%% check input
if(debugFlag == 1)
    tooth.vers = toothdata{1, 1};               % ��������
    tooth.tris = toothdata{2, 1};
    tooth.rVers = toothdata{9, 1};              % ������������
    tooth.rTris = toothdata{10, 1};
    writeOBJ('standardMesh.obj', aveBone.vers, aveBone.tris);
    objWriteVertices('dentalCenter.obj', dentalCenter);
    objWriteVertices('teethCenter.obj', glCenters);
    
    teeth.vers = [];
    teeth.tris = [];
    for i = 1:teethCount
        teeth.tris = [teeth.tris; toothdata{2, i} + size(teeth.vers,1)];
        teeth.vers = [teeth.vers; toothdata{1, i}];
    end
    writeOBJ('teeth.obj', teeth.vers, teeth.tris);
    allTeethCenter = mean(teeth.vers, 1);
    allGlCenter = mean(glCenters);
    objWriteVertices('allTeethCenter.obj', allTeethCenter);
    objWriteVertices('allGlCenter.obj', allGlCenter);
end
    

%% ����������ת�������xyƽ���ϲ��ñ�׼��Բ���

%1. �����任�����������ĵ�任���������ϵ��
teethCenterJC = bsxfun(@minus, glCenters, dentalCenter) * rotation;    

%2. �������ĵ������Բ������׼��Բ��û�������������ת��a*x^2 + c*y^2 + d*x + e*y + f = 0;
sampleVersJC = teethCenterJC(:, 1:2);
coff = fit_ellipse(sampleVersJC(:,1), sampleVersJC(:,2), 'standard');
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


% 3. �����Բ����ȡ����Ŀ��㣺����ѡȡ��������֮��ĵ�..
startVer = (sampleVersJC(1,:) + [-sampleVersJC(teethCount,1), sampleVersJC(teethCount,2)])/2;
idx1 = knnsearch(elliJC, startVer);                 % �����Բ�Ͼ���startVer����ĵ��������
idx2 = knnsearch(step, pi-step(idx1));
sel = round(linspace(idx1, idx2, 14));      % �����ھ��ȵ�ѡȡ14���㣻
targetVersJC = elliJC(sel,:);        
objWriteVertices('target1.obj', [elliJC(idx1, :), 0]);
objWriteVertices('target2.obj', [elliJC(idx2, :), 0]);

if(debugFlag == 1)
    objWriteVertices('sampleVers.obj', teethCenterJC);
    objWriteVertices('�����Բ.obj', [elliJC, zeros(size(elliJC, 1), 1)]);
    objWriteVertices('����Ŀ���.obj', [targetVersJC, zeros(size(targetVersJC, 1), 1)]);
end


%% ��׼�Ǳ��Ρ������Ƶ�λ�Ƶ�����Ŀ����λ�ã���Χ��������Ҳ����Ӧ��λ�ƣ�

% % ����Ȩ��
% conVersMat = reshape(conVers(:,1:2)', [1 2 14]); % ��ά����14�㣬ÿһ����һ�����Ƶ����ꣻ
% abVers2D = aveBone.vers(:,1:2);
% arrows = bsxfun(@minus, abVers2D, conVersMat);   % ��ά�������񶥵㵽ÿһ�����Ƶ��������
% arrowsLen = sum(arrows.^2, 2)/1000;
% weight = exp(-arrowsLen);
% % ��һ����
% weight = bsxfun(@rdivide, weight, sum(weight,3));    % rdivide(A, B) == A./B
% conArrows = (targetVersJC - conVers(:,1:2))';         % ���Ƶ�ָ����Բ������Ķ�ά���������洢�ɾ���
% conArrowsMat = reshape(conArrows, 1, 2, teethCount);
% offsetMat = bsxfun(@times, conArrowsMat, weight);
% % ���Ρ���ֻ��д��׼���۹Ƕ����xy����ֵ��
% aveBone.vers(:,1:2) = sum(offsetMat, 3) + aveBone.vers(:,1:2);
% % �任��global����ϵ��
% aveBone.vers = bsxfun(@plus, aveBone.vers * rotation', dentalCenter);


%% ���Բ�����ά�������Ȩ�أ�
versCount = size(aveBone.vers, 1);
arrowsLen = zeros(versCount, 14);
for i = 1: 14
    conVersMat = repmat(conVers(i, 1:2), versCount, 1);
    arrows = aveBone.vers(:, 1:2) - conVersMat;
    arrowsLen(:, i) = sum(arrows.^2, 2)/1000;
end

% ����ƫ��������Ȩ�أ�ÿ�����񶥵������ÿ�����Ƶ㶼��һ��Ȩ��ֵ���ܹ���versCount * 14����
weights = exp(-arrowsLen);

%   ��������һ����
sumCol = sum(weights, 2);
sumMat = repmat(sumCol, 1, 14);
weights = weights./sumMat;                  

conArrows = targetVersJC - conVers(:, 1:2);     % ���Ƶ�ָ�����Ŀ����������
offsetMatsArr = zeros(versCount, 2, 14);        % 14�㣬��i���ǵ�i�����Ƶ�ȷ�����������񶥵��ƫ��������
for i = 1:14
    arrow = conArrows(i, :);
    conArrowsMat = repmat(arrow, versCount, 1);
    offsetMatsArr(:, :, i) = bsxfun(@times, conArrowsMat, weights(:, i));
end
offsetArrows = sum(offsetMatsArr, 3);       % ���յ�ƫ������Ϊ14��ƫ�������ļӺͣ�
aveBone.vers(:,1:2) = aveBone.vers(:,1:2) + offsetArrows;

% �任��global����ϵ��
aveBone.vers = bsxfun(@plus, aveBone.vers * rotation', dentalCenter);


%% д������ݣ�
rootTeeth.vers = [];
rootTeeth.tris = [];
 
for i = 1:teethCount
    rootTeeth.tris = [rootTeeth.tris; toothdata{10,i} + size(rootTeeth.vers,1)];
    rootTeeth.vers = [rootTeeth.vers; toothdata{9,i}];
end


if(debugFlag == 1)
    writeOBJ('rootTeeth.obj', rootTeeth.vers, rootTeeth.tris);
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
 
 