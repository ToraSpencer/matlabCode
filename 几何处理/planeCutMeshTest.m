clc;
clear all;

% 平面切割网格，去除平面正侧的部分，保留负侧的部分；
vpDisThreshold = 1e-6;       % 点面距离小于此值时，认为点在面上；

[vers, tris] = readOBJfast('./data/tooth.obj'); 
versCount = size(vers, 1);
trisCount = size(tris, 1);
planeVert = [0, 0, 0];
planeNorm = [0, 0, 1];
    

%% 提取网格中在平面负侧的顶点：若(v - pv)*pn < vpDisThreshold则在负侧。
arrows = bsxfun(@minus, vers, planeVert);
vpDis = arrows * planeNorm';
tmpIdx = (1: versCount)';
negSideIdx = tmpIdx((vpDis < -vpDisThreshold), :);            % 平面负侧的顶点索引；
negSideVers = vers((vpDis < -vpDisThreshold), :);
inPlaneIdx = tmpIdx((abs(vpDis) <= vpDisThreshold), :);     % 平面上的顶点索引；
posiSideIdx =tmpIdx((vpDis > vpDisThreshold), :);            % 平面正侧的顶点索引；
objWriteVertices('negSideVers.obj', negSideVers);


%% 三角片矩阵每一个元素映射成对应顶点到平面的有向距离：
trisVec = reshape(tris', 3*trisCount, 1);       % 三角片数据拉成一条列向量；tris = (reshape(trisVec, 3, trisCount))';
triVers = vers(trisVec, :);
triVersArrows = bsxfun(@minus, triVers, planeVert);
triVersVPdis = triVersArrows * planeNorm';
triVersVPdisMat = (reshape(triVersVPdis, 3, trisCount))';

%% 若三角片中至少存在一个顶点在平面的负侧，则提取该三角片：
negFlagMat = (triVersVPdisMat < -vpDisThreshold);
flag1 = negFlagMat(:, 1) > 0;
flag2 = negFlagMat(:, 2) > 0;
flag3 = negFlagMat(:, 3) > 0;
flag = flag1 | flag2 | flag3;
tmpIdx = (1:trisCount)';
selectedTrisIdx = tmpIdx(flag);
[vers, tris] = selectedTris2Mesh(vers, tris, selectedTrisIdx);
versCount = size(vers, 1);
trisCount = size(tris, 1);
writeOBJ('selectedTrisMesh.obj', vers, tris);

%% shrinkEdges——找出所有穿刺边，沿着边方向位移不在平面负侧的那个顶点到平面上；
arrows = bsxfun(@minus, vers, planeVert);
vpDis = arrows * planeNorm';
tmpIdx = (1: versCount)';
naughtyVersIdx = tmpIdx((vpDis >= -vpDisThreshold), :);            % 不在平面负侧的顶点索引；

% 构建网格邻接关系
triEdges = [tris(:,2) tris(:,3); tris(:,3) tris(:,1); tris(:,1) tris(:,2)];    % 有向边，与三角片正法向对应；
triEdgesIdx = (1:3 * trisCount)';                   % 有向边的索引，列向量；
triEdgesCount = size(triEdgesIdx, 1);           % 有向边数；

adjSM_double = sparse(triEdges(:,1), triEdges(:,2), 1); % 顶点有向邻接矩阵（double)；下标(i, j)的元素为n == 存在n条有向边(i, j);          
adjSM = adjSM_double > 0;   % 顶点有向邻接矩阵（logical）；下标(i, j)的元素为1 == 存在有向边(i, j);

% 找出所有穿刺边：
flag1 = ismember(triEdges(:, 1), naughtyVersIdx);
flag2 = ismember(triEdges(:, 2), naughtyVersIdx);
flag = (flag1 &(~flag2)) | (~flag1 & flag2);
spearEdges = triEdges(flag, :);

% 穿刺边一律写成外点在前，内点在后的形式，然后按首顶点升序排列，去重；
spearEdgesInv = [spearEdges(:, 2), spearEdges(:, 1)];
flag = ismember(spearEdges(:, 2), naughtyVersIdx);              % 需要颠倒的边的索引；
spearEdges(flag, :) = spearEdgesInv(flag, :);
[~, tmp] = sort(spearEdges(:, 1));
tmpEdges = spearEdges(tmp, :);
spearEdges = unique(tmpEdges, 'rows');
objWriteEdges('spearEdges.obj', spearEdges, vers);

% for debug——测试planeEdgeIntersect()
e = spearEdges(1, :);
isctVer = planeEdgeIntersect(planeVert, planeNorm, e, vers);
objWriteVertices('eVers.obj', [isctVer; vers(e(1), :); vers(e(2), :)]);

disp('finished.');




