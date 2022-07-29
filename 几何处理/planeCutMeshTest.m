clc;
clear all;

% ƽ���и�����ȥ��ƽ������Ĳ��֣���������Ĳ��֣�
vpDisThreshold = 1e-6;       % �������С�ڴ�ֵʱ����Ϊ�������ϣ�

[vers, tris] = readOBJfast('./data/tooth.obj'); 
versCount = size(vers, 1);
trisCount = size(tris, 1);
planeVert = [0, 0, 0];
planeNorm = [0, 0, 1];
    

%% ��ȡ��������ƽ�渺��Ķ��㣺��(v - pv)*pn < vpDisThreshold���ڸ��ࡣ
arrows = bsxfun(@minus, vers, planeVert);
vpDis = arrows * planeNorm';
tmpIdx = (1: versCount)';
negSideIdx = tmpIdx((vpDis < -vpDisThreshold), :);            % ƽ�渺��Ķ���������
negSideVers = vers((vpDis < -vpDisThreshold), :);
inPlaneIdx = tmpIdx((abs(vpDis) <= vpDisThreshold), :);     % ƽ���ϵĶ���������
posiSideIdx =tmpIdx((vpDis > vpDisThreshold), :);            % ƽ������Ķ���������
objWriteVertices('negSideVers.obj', negSideVers);


%% ����Ƭ����ÿһ��Ԫ��ӳ��ɶ�Ӧ���㵽ƽ���������룺
trisVec = reshape(tris', 3*trisCount, 1);       % ����Ƭ��������һ����������tris = (reshape(trisVec, 3, trisCount))';
triVers = vers(trisVec, :);
triVersArrows = bsxfun(@minus, triVers, planeVert);
triVersVPdis = triVersArrows * planeNorm';
triVersVPdisMat = (reshape(triVersVPdis, 3, trisCount))';

%% ������Ƭ�����ٴ���һ��������ƽ��ĸ��࣬����ȡ������Ƭ��
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

%% shrinkEdges�����ҳ����д��̱ߣ����ű߷���λ�Ʋ���ƽ�渺����Ǹ����㵽ƽ���ϣ�
arrows = bsxfun(@minus, vers, planeVert);
vpDis = arrows * planeNorm';
tmpIdx = (1: versCount)';
naughtyVersIdx = tmpIdx((vpDis >= -vpDisThreshold), :);            % ����ƽ�渺��Ķ���������

% ���������ڽӹ�ϵ
triEdges = [tris(:,2) tris(:,3); tris(:,3) tris(:,1); tris(:,1) tris(:,2)];    % ����ߣ�������Ƭ�������Ӧ��
triEdgesIdx = (1:3 * trisCount)';                   % ����ߵ���������������
triEdgesCount = size(triEdgesIdx, 1);           % ���������

adjSM_double = sparse(triEdges(:,1), triEdges(:,2), 1); % ���������ڽӾ���double)���±�(i, j)��Ԫ��Ϊn == ����n�������(i, j);          
adjSM = adjSM_double > 0;   % ���������ڽӾ���logical�����±�(i, j)��Ԫ��Ϊ1 == ���������(i, j);

% �ҳ����д��̱ߣ�
flag1 = ismember(triEdges(:, 1), naughtyVersIdx);
flag2 = ismember(triEdges(:, 2), naughtyVersIdx);
flag = (flag1 &(~flag2)) | (~flag1 & flag2);
spearEdges = triEdges(flag, :);

% ���̱�һ��д�������ǰ���ڵ��ں����ʽ��Ȼ���׶����������У�ȥ�أ�
spearEdgesInv = [spearEdges(:, 2), spearEdges(:, 1)];
flag = ismember(spearEdges(:, 2), naughtyVersIdx);              % ��Ҫ�ߵ��ıߵ�������
spearEdges(flag, :) = spearEdgesInv(flag, :);
[~, tmp] = sort(spearEdges(:, 1));
tmpEdges = spearEdges(tmp, :);
spearEdges = unique(tmpEdges, 'rows');
objWriteEdges('spearEdges.obj', spearEdges, vers);

% for debug��������planeEdgeIntersect()
e = spearEdges(1, :);
isctVer = planeEdgeIntersect(planeVert, planeNorm, e, vers);
objWriteVertices('eVers.obj', [isctVer; vers(e(1), :); vers(e(2), :)]);

disp('finished.');




