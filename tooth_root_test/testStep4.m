clc;
clear all;

%% 4.连三角片 
   [patientCutVers, patientCutTris] = readOBJ('./step4/cutPatientTooth.obj');
   [rootCutVers, rootCutTris] = readOBJ('./step4/cutRootTooth.obj');
   axisPatient = readOBJ('./step4/axisPatient.obj');
    patientAxisTrans = axisPatient';
    
    
    % 4.1 确定牙冠、牙根边缘信息：
    hole = Calc_Boundary(patientCutTris);
    patientEdge = hole.boundary.edge;
    pEdgeVersIdx = patientEdge(:,1);
    pEdgeVers = patientCutVers(pEdgeVersIdx,:);
    
    hole = Calc_Boundary(rootCutTris);
    rootEdge = hole.boundary.edge;
    rEdgeVersIdx = rootEdge(:,1);
    rEdgeVers = rootCutVers(rEdgeVersIdx,:);
    
    % for debug
    OBJwriteVertices('pEdgeVers.obj', pEdgeVers);
    OBJwriteVertices('rEdgeVers.obj', rEdgeVers);
  
    % 4.2 边缘点集投影到二维平面，然后平滑
    mergeCenter = mean([pEdgeVers; rEdgeVers]);
    
    patientCircle = bsxfun(@minus, pEdgeVers, mergeCenter) * patientAxisTrans;
    patientCircle = smooth_loop(patientCircle(:,1:2), 0.01);      % 平滑操作？
 
 
    % 4.3 ？？？变形和位移？
    innerCircle = 0.5*bsxfun(@rdivide, patientCircle, normrow(patientCircle));

    rootCircle = bsxfun(@minus, rEdgeVers, mergeCenter) * patientAxisTrans;
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
    OBJwriteVertices('tempInner第一个点.obj', tempInner(1, :));
    OBJwriteVertices('tempOuter第一个点.obj', tempOuter(1, :));

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
        
    writeOBJ('融合网格test4.obj', finalVers, finalTris);
    
    
