 %% ����
 clc
 close all;
 clear all;
dbstop if error 


 %% 0 ׼������
    load('D:\workstation\gitRepositories\matlabCode\tooth_root_test\matlab.mat')
    
    order = 7;          % �������ݱ��
    
    toothdata = handle.model.toothdata;
    nT = size(toothdata,2);             % ������������
    
	% �������ĺ���������
    toothCenter = zeros(nT,3);
    for i = 1:nT
        toothCenter(i,:) = mean(toothdata{1,i});
    end
    dentalCenter = mean(toothCenter);
    
    if order >= 5 && order <=10        % ��Ӧ��1��2��3����
        [rootVers, rootTris] = readOBJ('root1.obj');
    else
        [rootVers, rootTris] = readOBJ('root2.obj');
    end
    patientVers = toothdata{1,order};
    patientTris = toothdata{2,order};
    patientDir = toothdata{7,order};

    writeOBJ('rootTooth.obj', rootVers, rootTris);
    writeOBJ('patientTooth.obj', patientVers, patientTris);

    % ���ݾֲ�����ϵ
    patientAxisTrans(:,3) = -patientDir';
    patientAxisTrans(:,1) = normalizerow(cross(toothCenter(order,:) - dentalCenter, -patientDir))';
    patientAxisTrans(:,2) = cross(patientAxisTrans(:,3)', patientAxisTrans(:,1));
    OBJwriteVertices('patientAxis.obj', patientAxisTrans');
    
    % ��ȡ��ƴ�����ݵ��������򡪡���ȡ���������������������ϵĲ��֡�
    patientCutVers = patientVers(1:toothdata{8,order},:);
    patientCutTris = patientTris(sum(patientTris <= toothdata{8,order}, 2)  == 3, :);

 
%% 1. �����и�����
    % ��ȡ��Ե���Ƶ�
    n = 20;
    hole = Calc_Boundary(patientCutTris);
    gumlineIdx = hole.boundary.edge(:,1);
    temp = ceil(linspace(1, length(gumlineIdx), n+1));
    pick = temp(1:end-1);           % ѡȡ��20�������ߵ㡣
    chosenVers = patientCutVers(gumlineIdx(pick),:);
    
    % for debug
    gumline = patientCutVers(gumlineIdx,:);
    OBJwriteVertices('gumline.obj', gumline);

    % ������������Ƶ�
    dis = 0.1;
    chosenVersMoved1 = bsxfun(@plus, chosenVers, dis*patientDir);
    chosenVersMoved2 = bsxfun(@minus, chosenVers, dis*patientDir);
    fittingPoints = [chosenVersMoved1; chosenVers; chosenVersMoved2];
    cv = [ones(n,1); zeros(n,1); -ones(n,1)];
 
    % ���������
    dist2 = pdist2(fittingPoints, fittingPoints);
    P = [ones(3*n,1), fittingPoints];
    A = [dist2, P; P', zeros(4,4)];
    b = [cv; zeros(4,1)];
    coeff = A\b;    % ���Է�����A*x == b�Ľ�������
    
    % ���ָ���
    syms x y z;
    f = coeff(1:3*n)' * ((x - fittingPoints(:,1)).^2 + (y - fittingPoints(:,2)).^2 + (z - fittingPoints(:,3)).^2).^(1/2) + ...
        [1 x y z] * coeff((3*n+1):end);
    f = matlabFunction(f);              % �ָ��溯����f(x,y,z) == 0;
    
    % for debug;
    OBJwriteVertices('fittingPoints.obj', fittingPoints);
 
    
%% 2. �и������
    patientValue = f(patientVers(:, 1), patientVers(:, 2), patientVers(:,3));
    cutPatientFlag = (patientValue < 0);
    patientCutVers2 = patientVers(cutPatientFlag, :);
    OBJwriteVertices('cutPatientVers2.obj', patientCutVers2);
    
    oldNewIdxInfo = -ones(size(patientVers, 1),1);
    reduceCount = 0;
    for i = 1:size(patientVers, 1)
        if(patientValue(i)>=0)
            reduceCount = reduceCount + 1;
            continue;
        end
        oldNewIdxInfo(i) = i-reduceCount;
    end
    patientTrisNewIdx = oldNewIdxInfo(patientTris);
    flag = (patientTrisNewIdx > 0);
    flag = sum(flag, 2);
    flag = (flag == 3);
    patientCutTris2 = patientTrisNewIdx(flag,:);
    writeOBJ('patientCut2.obj', patientCutVers2, patientCutTris2);
    patientCutVers = patientCutVers2;
    patientCutTris = patientCutTris2;           % �Լ��г����Ĳ������ڡ�

    
%% 3. ������Ե���ţ�������תƽ�ƣ��и�
    % ��׼��xy����
    Vctmp = bsxfun(@minus, patientCutVers, toothCenter(order,:)) * patientAxisTrans;
    max11 = max(Vctmp(:,1));
    min11 = min(Vctmp(:,1));
    max12 = max(rootVers(:,1)); 
    min12 = min(rootVers(:,1));
    scale1 = (max11 - min11)/(max12- min12);
    max21 = max(Vctmp(:,2));
    min21 = min(Vctmp(:,2));
    max22 = max(rootVers(:,2));
    min22 = min(rootVers(:,2));
    scale2 = (max21 - min21)/(max22 - min22);
    rootVers(:,1) = rootVers(:,1) * scale1;
    rootVers(:,2) = rootVers(:,2) * scale2;

    % ����׼��������ϵ����
    temp =  patientAxisTrans * rootVers';
    temp1 = rootVers * patientAxisTrans';
    rootVers = bsxfun(@plus, temp, toothCenter(order,:)')';

    rootVers2 = temp' + repmat(toothCenter(order, :), size(temp', 1),1);
    rootVers3 = temp1 + repmat(toothCenter(order, :), size(temp1, 1),1);
    
    % ��ȡ��׼������������
    isRoot = f(rootVers(:,1),rootVers(:,2),rootVers(:,3)) > 0;
    Fp = tt(rootTris);
    tag = double( sum(isRoot(rootTris),2) == 3 );
    C = connected_region(Fp, tag);      % ĳһ������Ƭ��ɢ�γ���ͨͼ�������ж��ͼ��ȡ�����Ǹ���
    [~, ind] = sort( cellfun('length', C), 'descend' );
    [rootCutVers, rootCutTris] = Remove_Point(rootVers, rootTris(C{ind(1)},:));

    writeOBJ('cutPatientTooth.obj', patientCutVers, patientCutTris);
    writeOBJ('cutRootTooth.obj', rootCutVers, rootCutTris);

    
%% 4.������Ƭ 

    % ��ƽ���Ͻ������ӹ�ϵ
    hole = Calc_Boundary(patientCutTris);
    patientEdge = hole.boundary.edge;
    pEdgeVersIdx = patientEdge(:,1);
    pEdgeVers = patientCutVers(pEdgeVersIdx,:);
    
    hole = Calc_Boundary(rootCutTris);
    rootEdge = hole.boundary.edge;
    rEdgeVersIdx = rootEdge(:,1);
    rEdgeVers = rootCutVers(rEdgeVersIdx,:);
  
    mergeCenter = mean([pEdgeVers; rEdgeVers]);
    
    innerCircle = bsxfun(@minus, pEdgeVers, mergeCenter) * patientAxisTrans;
    innerCircle = smooth_loop(innerCircle(:,1:2), 0.01);      % ƽ��������
    innerCircle = 0.5*bsxfun(@rdivide, innerCircle, normrow(innerCircle));
    
    outerCircle = bsxfun(@minus, rEdgeVers, mergeCenter) * patientAxisTrans;
    outerCircle = smooth_loop(outerCircle(:,1:2), 0.01);      
    outerCircle = bsxfun(@rdivide, outerCircle, normrow(outerCircle));

    innerCount = size(innerCircle,1);
    outerCount = size(outerCircle,1);
    edgeInPlane = [[1:innerCount; [2:innerCount,1]]'; [1:outerCount; [2:outerCount,1]]' + innerCount];
    versInPlane = [innerCircle; outerCircle];
    [~, trisInPlane] = triangle(versInPlane, edgeInPlane, mean(innerCircle), 'NoBoundarySteiners');
    %Fa = Fa(:, [3,2,1]);
    
    % ���ں������ں�
    finalVers = [patientCutVers; rootCutVers];
    edgeVersIdx = [pEdgeVersIdx; rEdgeVersIdx + size(patientCutVers,1)];
    addTris = edgeVersIdx(trisInPlane);       % ��ά�ռ������������ɵ�Ƭת��Ϊ�ں������е�����Ƭ��
    finalTris = [patientCutTris; rootCutTris + size(patientCutVers,1); addTris];

    
%% 5.������������״��
    % �ںϲ��ֹ⻬
    A = adjacency_matrix(finalTris);
    L = A - diag(sparse(sum(A,2)));
    lhs = L;
    L2 = L*L;
    lhs(rEdgeVersIdx + size(patientCutVers,1), :) = L2(rEdgeVersIdx + size(patientCutVers,1), :);
    rhs = L*finalVers;
    rhs(rEdgeVersIdx + size(patientCutVers,1), :) = 0;
    finalVers = solve_equation(lhs, rhs, 1:size(patientCutVers,1), patientCutVers);
 
    
    writeOBJ('���ս��.obj', finalVers, finalTris);

    disp('finished');
    
%%
    figure
    drawMesh(versInPlane, trisInPlane);

