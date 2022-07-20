clc;
clear all;
close all;
debugFlag = 1;

%% 读输入数据，prepare data:
load('inputData/data.mat');

aveBone.vers = [];
aveBone.tris = [];

[aveBone.vers, aveBone.tris] = readOBJ('inputData/aveBone_lower_stan.obj');       % 标准下颌牙槽骨网格；
conVers = [                 % 标准牙槽骨的控制点。
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
rotation = handle.dentalFrame;
toothdata = handle.model.toothdata;

teethCount = size(toothdata,2);

% 牙齿中心
teethCenter = zeros(teethCount,3);
for i = 1:teethCount
    teethCenter(i,:) = mean(toothdata{1,i}(toothdata{3,i},:));          % 字段三是牙龈线点在牙齿网格中的索引；
end


%% check input
if(debugFlag == 1)
    objWriteVertices('conVers.obj', conVers);
    tooth.vers = toothdata{1, 1};               % 牙齿网格
    tooth.tris = toothdata{2, 1};
    tooth.rVers = toothdata{9, 1};              % 带根牙齿网格；
    tooth.rTris = toothdata{10, 1};
    writeOBJ('tooth1.obj', tooth.vers, tooth.tris);
    writeOBJ('rootTooth1.obj', tooth.rVers, tooth.rTris);
    objWriteVertices('dentalCenter.obj', dentalCenter);
    objWriteVertices('teethCenter.obj', teethCenter);
end
    

%% 将牙齿中心转换到牙颌xy平面上并用标准椭圆拟合

%1. 逆仿射变换，将牙齿中心点变换到牙颌坐标系；
teethCenterJC = bsxfun(@minus, teethCenter, dentalCenter) * rotation;    

%2. 牙齿中心点拟合椭圆——标准椭圆（没有相对坐标轴旋转）a*x^2 + c*y^2 + d*x + e*y + f = 0;
sampleVersJC = teethCenterJC(:, 1:2);
coff = fit_ellipse(sampleVersJC(:,1), sampleVersJC(:,2), 'standard');
a = coff(1);
c = coff(2);
d = coff(3);
e = coff(4);
f = coff(5);

center = [-d/(2*a), -e/(2*c)];            % 椭圆中心坐标；
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


% 3. 拟合椭圆中提取控制目标点：均匀选取两端牙齿之间的点..
startVer = (sampleVersJC(1,:) + [-sampleVersJC(teethCount,1), sampleVersJC(teethCount,2)])/2;
idx1 = knnsearch(elliJC, startVer);                 % 拟合椭圆上距离startVer最近的点的索引；
idx2 = knnsearch(step, pi-step(idx1));
sel = round(linspace(idx1, idx2, 14));      % 区间内均匀地选取14个点；
targetVersJC = elliJC(sel,:);        
objWriteVertices('target1.obj', [elliJC(idx1, :), 0]);
objWriteVertices('target2.obj', [elliJC(idx2, :), 0]);

if(debugFlag == 1)
    objWriteVertices('sampleVers.obj', teethCenterJC);
    objWriteVertices('拟合椭圆.obj', [elliJC, zeros(size(elliJC, 1), 1)]);
    objWriteVertices('控制目标点.obj', [targetVersJC, zeros(size(targetVersJC, 1), 1)]);
end


%% 标准颌骨变形——控制点位移到控制目标点的位置，周围其他顶点也做相应的位移；

% % 计算权重
% conVersMat = reshape(conVers(:,1:2)', [1 2 14]); % 三维矩阵，14层，每一层是一个控制点坐标；
% arrows = bsxfun(@minus, aveBone.vers(:,1:2), conVersMat);   % 三维矩阵，网格顶点到每一个控制点的向量；
% arrowsLen = sum(arrows.^2, 2)/1000;
% weight = exp(-arrowsLen);
% % 归一化：
% weight = bsxfun(@rdivide, weight, sum(weight,3));    % rdivide(A, B) == A./B
% conArrows = (targetVersJC - conVers(:,1:2))';         % 控制点指向椭圆采样点的二维列向量，存储成矩阵；
% conArrowsMat = reshape(conArrows, 1, 2, teethCount);
% offsetMat = bsxfun(@times, conArrowsMat, weight);
% % 变形——只改写标准牙槽骨顶点的xy坐标值；
% aveBone.vers(:,1:2) = sum(offsetMat, 3) + aveBone.vers(:,1:2);
% % 变换到global坐标系；
% aveBone.vers = bsxfun(@plus, aveBone.vers * dentalFrame', dentalCenter);


%% 尝试不用三维矩阵计算权重：
versCount = size(aveBone.vers, 1);
arrowsLen = zeros(versCount, 14);
for i = 1: 14
    conVersMat = repmat(conVers(i, 1:2), versCount, 1);
    arrows = aveBone.vers(:, 1:2) - conVersMat;
    arrowsLen(:, i) = sum(arrows.^2, 2)/1000;
end

% 计算偏移向量的权重，每个网格顶点相对于每个控制点都有一个权重值，总共有versCount * 14个；
weights = exp(-arrowsLen);

%   列向量归一化：
sumCol = sum(weights, 2);
sumMat = repmat(sumCol, 1, 14);
weights = weights./sumMat;                  

conArrows = targetVersJC - conVers(:, 1:2);     % 控制点指向控制目标点的向量；
offsetMatsArr = zeros(versCount, 2, 14);        % 14层，第i层是第i个控制点确定的所有网格顶点的偏移向量；
for i = 1:14
    arrow = conArrows(i, :);
    conArrowsMat = repmat(arrow, versCount, 1);
    offsetMatsArr(:, :, i) = bsxfun(@times, conArrowsMat, weights(:, i));
end
offsetArrows = sum(offsetMatsArr, 3);       % 最终的偏移向量为14个偏移向量的加和；
aveBone.vers(:,1:2) = aveBone.vers(:,1:2) + offsetArrows;

% 变换到global坐标系；
aveBone.vers = bsxfun(@plus, aveBone.vers * rotation', dentalCenter);


%% 写输出数据：
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
 
 