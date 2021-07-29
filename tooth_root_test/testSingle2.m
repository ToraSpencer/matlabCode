 clc
 close all;
 clear all;
dbstop if error 


 %% 0 准备数据
[patientCutVers, patientCutTris] = readOBJ('切割病人牙冠.obj');
[rootCutVers, rootCutTris] = readOBJ('切割标准牙根.obj');
patientAxis = readOBJ('axisPatient.obj');
patientCenter = mean(patientCutVers, 1);

 
zdir = patientAxis(2, :);
ydir = patientAxis(3, :);
xdir = cross(ydir, zdir);


printDir('xdir.obj', xdir, patientCenter);
printDir('ydir.obj', ydir, patientCenter);
printDir('zdir.obj', zdir, patientCenter);
patientAxisTrans = [xdir; ydir; zdir]';


    
%% 4.连三角片 

    % 4.1 确定牙冠、牙根边缘信息：
    hole = Calc_Boundary(patientCutTris);
    patientEdge = hole.boundary.edge;
    pEdgeVersIdx = patientEdge(:,1);
    pEdgeVers = patientCutVers(pEdgeVersIdx,:);
    
    hole = Calc_Boundary(rootCutTris);
    rootEdge = hole.boundary.edge;
    rEdgeVersIdx = rootEdge(:,1);
    rEdgeVers = rootCutVers(rEdgeVersIdx,:);
    
    %   for debug:
    OBJwriteVertices('rEdgeVers.obj', rEdgeVers);
    OBJwriteVertices('pEdgeVers.obj', pEdgeVers);
  
    % 4.2 边缘点集投影到二维平面，然后平滑
    mergeCenter = mean([pEdgeVers; rEdgeVers]);
    
    patientCircle = bsxfun(@minus, pEdgeVers, mergeCenter) * patientAxisTrans;
    patientCircle = smooth_loop(patientCircle(:,1:2), 0.001);      % 平滑操作？
    rootCircle = bsxfun(@minus, rEdgeVers, mergeCenter) * patientAxisTrans;
    

    % 4.3 ？？？变形和位移？
    innerCircle = 0.5*bsxfun(@rdivide, patientCircle, normrow(patientCircle));


    rootCircle = smooth_loop(rootCircle(:,1:2), 0.01);    
    temp = normrow(rootCircle);
    temp = repmat(temp,1, 2);
    temp2 = rootCircle./temp;
    outerCircle = bsxfun(@rdivide, rootCircle, normrow(rootCircle));
    
    % for debug
    tempOuter = [outerCircle, zeros(size(outerCircle, 1), 1)];
    tempInner = [innerCircle, zeros(size(innerCircle, 1), 1)];
    OBJwriteVertices('tempInner.obj', tempInner);
    OBJwriteVertices('tempOuter.obj', tempOuter);
     OBJwriteVertices('twoCircles.obj', [tempOuter; tempInner]);

    % 4.4 二维点集三角剖分
    innerCount = size(innerCircle,1);
    outerCount = size(outerCircle,1);
    edgeInPlane = [[1:innerCount; [2:innerCount,1]]'; [1:outerCount; [2:outerCount,1]]' + innerCount];
    versInPlane = [innerCircle; outerCircle];
    [~, trisInPlane] = triangle(versInPlane, edgeInPlane, mean(innerCircle), 'NoBoundarySteiners');
    figure
    drawMesh(versInPlane, trisInPlane);
    
    % 4.5 牙冠和牙根融合
    finalVers = [patientCutVers; rootCutVers];
    rEdgeVersIdx = rEdgeVersIdx + size(patientCutVers,1);
    edgeVersIdx = [pEdgeVersIdx; rEdgeVersIdx];
    addTris = edgeVersIdx(trisInPlane);       % 二维空间中两个环连成的片转化为融合网格中的三角片。
    finalTris = [patientCutTris; rootCutTris + size(patientCutVers,1); addTris];
        
    writeOBJ('融合网格.obj', finalVers, finalTris);
    
    
%% 5.修整牙根的形状：
    % 融合部分光滑
    adjMat = adjacency_matrix(finalTris);   % 融合网格的邻接矩阵,维度为versCount*versCount，两个顶点有边连接则对应元素为1；
 
    temp = sum(adjMat, 2);   
    L = adjMat - diag(sparse(sum(adjMat,2)));
 
    %for debug
    flag = (temp == zeros(size(temp, 2), 1));
    
    A = L;
    
    L2 = L*L;
    A(rEdgeVersIdx, :) = L2(rEdgeVersIdx, :);
 
    B = L*finalVers;
    B(rEdgeVersIdx, :) = 0;
    
    Acon = 1:size(patientCutVers,1);
    Bcon = patientCutVers;
    finalVers = solve_equation_modified(A, B, Acon, Bcon);
    
    writeOBJ('最终结果.obj', finalVers, finalTris);

    disp('finished');
    
    
%%


