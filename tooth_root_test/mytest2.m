 %% 单个
 clc
 close all;
 clear all;
dbstop if error 


 %% 0 准备数据
    load('D:\workstation\gitRepositories\matlabCode\tooth_root_test\matlab.mat')
    
    order = 7;          % 病人牙齿编号
    
    toothdata = handle.model.toothdata;
    nT = size(toothdata,2);             % 病人牙齿总数
    
	% 牙齿中心和牙合中心
    toothCenter = zeros(nT,3);
    for i = 1:nT
        toothCenter(i,:) = mean(toothdata{1,i});
    end
    dentalCenter = mean(toothCenter);
    
    if order >= 5 && order <=10        % 对应着1，2，3号牙
        [rootVers, rootTris] = readOBJ('root1.obj');
    else
        [rootVers, rootTris] = readOBJ('root2.obj');
    end
    patientVers = toothdata{1,order};
    patientTris = toothdata{2,order};
    patientDir = toothdata{7,order};

    writeOBJ('rootTooth.obj', rootVers, rootTris);
    writeOBJ('patientTooth.obj', patientVers, patientTris);

    % 牙齿局部坐标系
    patientAxisTrans(:,3) = -patientDir';
    patientAxisTrans(:,1) = normalizerow(cross(toothCenter(order,:) - dentalCenter, -patientDir))';
    patientAxisTrans(:,2) = cross(patientAxisTrans(:,3)', patientAxisTrans(:,1));
    OBJwriteVertices('patientAxis.obj', patientAxisTrans');
    
    % 提取待拼接牙齿的牙冠区域——提取病人牙齿网格牙龈线以上的部分。
    patientCutVers = patientVers(1:toothdata{8,order},:);
    patientCutTris = patientTris(sum(patientTris <= toothdata{8,order}, 2)  == 3, :);

 
%% 1. 计算切割曲面
    % 提取边缘控制点
    n = 20;
    hole = Calc_Boundary(patientCutTris);
    gumlineIdx = hole.boundary.edge(:,1);
    temp = ceil(linspace(1, length(gumlineIdx), n+1));
    pick = temp(1:end-1);           % 选取的20个牙龈线点。
    chosenVers = patientCutVers(gumlineIdx(pick),:);
    
    % for debug
    gumline = patientCutVers(gumlineIdx,:);
    OBJwriteVertices('gumline.obj', gumline);

    % 径向基函数控制点
    dis = 0.1;
    chosenVersMoved1 = bsxfun(@plus, chosenVers, dis*patientDir);
    chosenVersMoved2 = bsxfun(@minus, chosenVers, dis*patientDir);
    fittingPoints = [chosenVersMoved1; chosenVers; chosenVersMoved2];
    cv = [ones(n,1); zeros(n,1); -ones(n,1)];
 
    % 径向基函数
    dist2 = pdist2(fittingPoints, fittingPoints);
    P = [ones(3*n,1), fittingPoints];
    A = [dist2, P; P', zeros(4,4)];
    b = [cv; zeros(4,1)];
    coeff = A\b;    % 线性方程组A*x == b的解向量；
    
    % 画分割面
    syms x y z;
    f = coeff(1:3*n)' * ((x - fittingPoints(:,1)).^2 + (y - fittingPoints(:,2)).^2 + (z - fittingPoints(:,3)).^2).^(1/2) + ...
        [1 x y z] * coeff((3*n+1):end);
    f = matlabFunction(f);              % 分割面函数：f(x,y,z) == 0;
    
    % for debug;
    OBJwriteVertices('fittingPoints.obj', fittingPoints);
 
    
%% 2. 切割病人牙齿
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
    patientCutTris = patientCutTris2;           % 自己切出来的病人牙冠。

    
%% 3. 牙根边缘缩放，整体旋转平移，切割
    % 标准牙xy缩放
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

    % 将标准牙的坐标系对齐
    temp =  patientAxisTrans * rootVers';
    temp1 = rootVers * patientAxisTrans';
    rootVers = bsxfun(@plus, temp, toothCenter(order,:)')';

    rootVers2 = temp' + repmat(toothCenter(order, :), size(temp', 1),1);
    rootVers3 = temp1 + repmat(toothCenter(order, :), size(temp1, 1),1);
    
    % 提取标准牙的牙根区域
    isRoot = f(rootVers(:,1),rootVers(:,2),rootVers(:,3)) > 0;
    Fp = tt(rootTris);
    tag = double( sum(isRoot(rootTris),2) == 3 );
    C = connected_region(Fp, tag);      % 某一个三角片扩散形成连通图，可能有多个图，取最大的那个。
    [~, ind] = sort( cellfun('length', C), 'descend' );
    [rootCutVers, rootCutTris] = Remove_Point(rootVers, rootTris(C{ind(1)},:));

    writeOBJ('cutPatientTooth.obj', patientCutVers, patientCutTris);
    writeOBJ('cutRootTooth.obj', rootCutVers, rootCutTris);

    
%% 4.连三角片 

    % 在平面上建立连接关系
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
    innerCircle = smooth_loop(innerCircle(:,1:2), 0.01);      % 平滑操作？
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
    
    % 牙冠和牙根融合
    finalVers = [patientCutVers; rootCutVers];
    edgeVersIdx = [pEdgeVersIdx; rEdgeVersIdx + size(patientCutVers,1)];
    addTris = edgeVersIdx(trisInPlane);       % 二维空间中两个环连成的片转化为融合网格中的三角片。
    finalTris = [patientCutTris; rootCutTris + size(patientCutVers,1); addTris];

    
%% 5.修整牙根的形状：
    % 融合部分光滑
    A = adjacency_matrix(finalTris);
    L = A - diag(sparse(sum(A,2)));
    lhs = L;
    L2 = L*L;
    lhs(rEdgeVersIdx + size(patientCutVers,1), :) = L2(rEdgeVersIdx + size(patientCutVers,1), :);
    rhs = L*finalVers;
    rhs(rEdgeVersIdx + size(patientCutVers,1), :) = 0;
    finalVers = solve_equation(lhs, rhs, 1:size(patientCutVers,1), patientCutVers);
 
    
    writeOBJ('最终结果.obj', finalVers, finalTris);

    disp('finished');
    
%%
    figure
    drawMesh(versInPlane, trisInPlane);

