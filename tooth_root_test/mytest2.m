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

    
%%  4.������Ƭ 

    % ��ƽ���Ͻ������ӹ�ϵ
    hole = Calc_Boundary(patientCutTris);
    idx_i = hole.boundary.edge(:,1);
    p_i_3d = patientCutVers(idx_i,:);
    hole = Calc_Boundary(rootCutTris);
    idx_o = hole.boundary.edge(:,1);
    p_o_3d = rootCutVers(idx_o,:);
    proj = zeros(3,3);
    proj(:,3) = patientDir';
    proj(:,1) = normalizerow(cross(proj(:,3)', [0 0 1]))';
    proj(:,2) = cross(proj(:,3)', proj(:,1)')';
    p_i_2d = bsxfun(@minus, p_i_3d, mean([p_i_3d; p_o_3d])) * proj;
    p_i_2d = smooth_loop(p_i_2d(:,1:2), 0.01);
    p_i_2d = 0.5*bsxfun(@rdivide, p_i_2d, normrow(p_i_2d));
    p_o_2d = bsxfun(@minus, p_o_3d, mean([p_i_3d; p_o_3d])) * proj;
    p_o_2d = smooth_loop(p_o_2d(:,1:2), 0.01);
    p_o_2d = bsxfun(@rdivide, p_o_2d, normrow(p_o_2d));

    n_i = size(p_i_2d,1);
    n_o = size(p_o_2d,1);
    E = [[1:n_i; [2:n_i,1]]'; [1:n_o; [2:n_o,1]]' + n_i];
    [~, Fa] = triangle([p_i_2d; p_o_2d], E, mean(p_i_2d), 'NoBoundarySteiners');
    Fa = Fa(:, [3,2,1]);

%% 5.������������״��

    % ���ں������ں�
    finalVers = [patientCutVers; rootCutVers];
    reidx = [idx_i; idx_o + size(patientCutVers,1)];
    finalTris = [patientCutTris; rootCutTris + size(patientCutVers,1); reidx(Fa)];

    % �ںϲ��ֹ⻬
    A = adjacency_matrix(finalTris);
    L = A - diag(sparse(sum(A,2)));
    lhs = L;
    L2 = L*L;
    lhs(idx_o + size(patientCutVers,1), :) = L2(idx_o + size(patientCutVers,1), :);
    rhs = L*finalVers;
    rhs(idx_o + size(patientCutVers,1), :) = 0;
    finalVers = solve_equation(lhs, rhs, 1:size(patientCutVers,1), patientCutVers);
 
    
    writeOBJ('���ս��.obj', finalVers, finalTris);

    disp('finished');

