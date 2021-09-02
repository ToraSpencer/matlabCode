 %% 单个
 clc
 close all;
 clear all;
dbstop if error 


 %% 0 准备数据
[patientCutVers, patientCutTris] = readOBJ('切割病人牙冠.obj');
[rootCutVers, rootCutTris] = readOBJ('切割标准牙根.obj');
patientAxis = readOBJ('axisPatient.obj');
patientCenter = mean(patientCutVers, 1);

% 		std::swap(axisPatient[1], axisPatient[2]);
% 		axisPatient[0] = axisPatient[1].Cross(axisPatient[2]);
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
  
    % 4.2 边缘点集投影到二维平面
    mergeCenter = mean([pEdgeVers; rEdgeVers]);
    
    patientCircle = bsxfun(@minus, pEdgeVers, mergeCenter) * patientAxisTrans;
    rootCircle = bsxfun(@minus, rEdgeVers, mergeCenter) * patientAxisTrans;
    
    OBJwriteVertices('patientCircle.obj', patientCircle);
    patientCircleProj = [patientCircle(:, 1:2),  zeros(size(patientCircle, 1), 1)];
    OBJwriteVertices('patientCircleProj.obj', patientCircleProj);
    
    innerCircle = patientCircle(:, 1:2);
    innerCount = size(innerCircle, 1);
    outerCircle = rootCircle(:, 1:2);
    outerCount = size(outerCircle, 1);
  
    % 转换到极坐标系：
    innerCirclePolar = [normrow(innerCircle), acos(innerCircle(:,1)./normrow(innerCircle))];
    for i = 1:innerCount
        if innerCircle(i, 2) < 0
            innerCirclePolar(i, 2) = -innerCirclePolar(i, 2);
        end
    end
    theta1 = innerCirclePolar(1, 2);
    theta2 = innerCirclePolar(round(innerCount/4), 2);
    temp = handleArc(theta2 - theta1); 
    if (temp >= 0 && temp < pi)
        deltaTheta = 2*pi/innerCount;
    else
        deltaTheta = -2*pi/innerCount;
    end
    
    for i = 1:innerCount
        innerCircle(i, 1) = 0.5*cos(theta1 + deltaTheta*(i-1));
        innerCircle(i, 2) = 0.5*sin(theta1 + deltaTheta*(i-1));
    end
    
 
    outerCirclePolar = [normrow(outerCircle), acos(outerCircle(:,1)./normrow(outerCircle))];
    for i = 1:outerCount
        if outerCircle(i, 2) < 0
            outerCirclePolar(i, 2) = -outerCirclePolar(i, 2);
        end
    end
    theta1 = outerCirclePolar(1, 2);
    theta2 = outerCirclePolar(round(outerCount/4), 2);
    temp = handleArc(theta2 - theta1); 
    if (temp >= 0 && temp < pi)
        deltaTheta = 2*pi/outerCount;
    else
        deltaTheta = -2*pi/outerCount;
    end
    
    for i = 1:outerCount
        outerCircle(i, 1) = 1.0*cos(theta1 + deltaTheta*(i-1));
        outerCircle(i, 2) = 1.0*sin(theta1 + deltaTheta*(i-1));
    end
 
    % for debug
    tempOuter = [outerCircle, zeros(size(outerCircle, 1), 1)];
    tempInner = [innerCircle, zeros(size(innerCircle, 1), 1)];
    OBJwriteVertices('映射为标准圆的tempInner.obj', tempInner);
    OBJwriteVertices('映射为标准圆的tempOuter.obj', tempOuter);
    OBJwriteVertices('twoCircles.obj', [tempOuter; tempInner]);

    % 4.4 二维点集三角剖分
    outerCount = size(outerCircle,1);
    edgeInPlane = [[1:innerCount; [2:innerCount,1]]'; [1:outerCount; [2:outerCount,1]]' + innerCount];
    versInPlane = [innerCircle; outerCircle];
    [~, trisInPlane] = triangle(versInPlane, edgeInPlane, mean(innerCircle), 'NoBoundarySteiners');

    % for debug
    versInPlaneTemp = zeros(size(versInPlane, 1), 3);
    versInPlaneTemp(:, 1 : 2) = versInPlane;
    writeOBJ('三角剖分生成的圆环.obj', versInPlaneTemp, trisInPlane);
    
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

