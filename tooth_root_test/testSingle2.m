 clc
 close all;
 clear all;
dbstop if error 


 %% 0 ׼������
[patientCutVers, patientCutTris] = readOBJ('�и������.obj');
[rootCutVers, rootCutTris] = readOBJ('�и��׼����.obj');
patientAxis = readOBJ('axisPatient.obj');
patientCenter = mean(patientCutVers, 1);

 
zdir = patientAxis(2, :);
ydir = patientAxis(3, :);
xdir = cross(ydir, zdir);


printDir('xdir.obj', xdir, patientCenter);
printDir('ydir.obj', ydir, patientCenter);
printDir('zdir.obj', zdir, patientCenter);
patientAxisTrans = [xdir; ydir; zdir]';


    
%% 4.������Ƭ 

    % 4.1 ȷ�����ڡ�������Ե��Ϣ��
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
  
    % 4.2 ��Ե�㼯ͶӰ����άƽ�棬Ȼ��ƽ��
    mergeCenter = mean([pEdgeVers; rEdgeVers]);
    
    patientCircle = bsxfun(@minus, pEdgeVers, mergeCenter) * patientAxisTrans;
    patientCircle = smooth_loop(patientCircle(:,1:2), 0.001);      % ƽ��������
    rootCircle = bsxfun(@minus, rEdgeVers, mergeCenter) * patientAxisTrans;
    

    % 4.3 ���������κ�λ�ƣ�
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

    % 4.4 ��ά�㼯�����ʷ�
    innerCount = size(innerCircle,1);
    outerCount = size(outerCircle,1);
    edgeInPlane = [[1:innerCount; [2:innerCount,1]]'; [1:outerCount; [2:outerCount,1]]' + innerCount];
    versInPlane = [innerCircle; outerCircle];
    [~, trisInPlane] = triangle(versInPlane, edgeInPlane, mean(innerCircle), 'NoBoundarySteiners');
    figure
    drawMesh(versInPlane, trisInPlane);
    
    % 4.5 ���ں������ں�
    finalVers = [patientCutVers; rootCutVers];
    rEdgeVersIdx = rEdgeVersIdx + size(patientCutVers,1);
    edgeVersIdx = [pEdgeVersIdx; rEdgeVersIdx];
    addTris = edgeVersIdx(trisInPlane);       % ��ά�ռ������������ɵ�Ƭת��Ϊ�ں������е�����Ƭ��
    finalTris = [patientCutTris; rootCutTris + size(patientCutVers,1); addTris];
        
    writeOBJ('�ں�����.obj', finalVers, finalTris);
    
    
%% 5.������������״��
    % �ںϲ��ֹ⻬
    adjMat = adjacency_matrix(finalTris);   % �ں�������ڽӾ���,ά��ΪversCount*versCount�����������б��������ӦԪ��Ϊ1��
 
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
    
    writeOBJ('���ս��.obj', finalVers, finalTris);

    disp('finished');
    
    
%%


