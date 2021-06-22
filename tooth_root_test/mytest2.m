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
        nameStr = ['patientAxis', num2str(i), '.obj'];
        currentAxis = toothdata{7,i};
        OBJwriteVertices(nameStr, currentAxis);
    end
    dentalCenter = mean(toothCenter);
    
    if order >= 5 && order <=10        % 对应着1，2，3号牙
        [rootVers, rootTris] = readOBJ('root1.obj');
    else
        [rootVers, rootTris] = readOBJ('root2.obj');
    end
    patientVers = toothdata{1,order};
    patientTris = toothdata{2,order};
    patientAxis = toothdata{7,order};

    writeOBJ('rootTooth.obj', rootVers, rootTris);
    writeOBJ('patientTooth.obj', patientVers, patientTris);

    % 牙齿局部坐标系
    toothCoord(:,1) = normalizerow(cross(toothCenter(order,:) - dentalCenter, -patientAxis))';
    toothCoord(:,3) = -patientAxis';
    toothCoord(:,2) = cross(toothCoord(:,3)', toothCoord(:,1));

    
    % 提取待拼接牙齿的牙冠区域——提取病人牙齿网格牙龈线以上的部分。
    patientCutVers = patientVers(1:toothdata{8,order},:);
    patientCutTris = patientTris(sum(patientTris <= toothdata{8,order}, 2)  == 3, :);

 
    
%% 1. 计算切割曲面
    % 提取边缘控制点
    n = 20;
    hole = Calc_Boundary(patientCutTris);
    edgeVersIdx = hole.boundary.edge(:,1);
    temp = ceil(linspace(1, length(edgeVersIdx), n+1));
    pick = temp(1:end-1);
    cp0 = patientCutVers(edgeVersIdx(pick),:);
    
    % for debug
    gumline = patientCutVers(edgeVersIdx,:);
    OBJwriteVertices('gumline.obj', gumline);

    % 径向基函数控制点
    eta = 0.1;
    cp = [bsxfun(@plus, cp0, eta*patientAxis); cp0; bsxfun(@minus, cp0, eta*patientAxis)];
    cv = [ones(n,1); zeros(n,1); -ones(n,1)];

    % 径向基函数
    dist2 = pdist2(cp, cp);
    P = [ones(3*n,1), cp];
    A = [dist2, P; P', zeros(4,4)];
    B = [cv; zeros(4,1)];
    coeff = A\B;

%% 2.牙根切割边缘缩放，整体旋转平移

    % 标准牙xy缩放
    Vctmp = bsxfun(@minus, patientCutVers, toothCenter(order,:)) * toothCoord;
    rootVers(:,1) = rootVers(:,1) * (max(Vctmp(:,1)) - min(Vctmp(:,1)))/(max(rootVers(:,1)) - min(rootVers(:,1)));
    rootVers(:,2) = rootVers(:,2) * (max(Vctmp(:,2)) - min(Vctmp(:,2)))/(max(rootVers(:,2)) - min(rootVers(:,2)));

    % 将标准牙的坐标系对齐
    rootVers = bsxfun(@plus, toothCoord * rootVers', toothCenter(order,:)')';

    % 画分割面
    syms x y z;
    f = coeff(1:3*n)' * ((x - cp(:,1)).^2 + (y - cp(:,2)).^2 + (z - cp(:,3)).^2).^(1/2) + ...
        [1 x y z] * coeff((3*n+1):end);
    f = matlabFunction(f);              % 分割面函数：f(x,y,z) == 0;


    % 提取标准牙的牙根区域
    isRoot = f(rootVers(:,1),rootVers(:,2),rootVers(:,3)) > 0;
    Fp = tt(rootTris);
    tag = double( sum(isRoot(rootTris),2) == 3 );
    C = connected_region(Fp, tag);
    [~, ind] = sort( cellfun('length', C), 'descend' );
    [rootCutVers, rootCutTris] = Remove_Point(rootVers, rootTris(C{ind(1)},:));

    writeOBJ('cutPatientTooth.obj', patientCutVers, patientCutTris);
    writeOBJ('cutRootTooth.obj', rootCutVers, rootCutTris);

    
%%  3.连三角片 

    % 在平面上建立连接关系
    hole = Calc_Boundary(patientCutTris);
    idx_i = hole.boundary.edge(:,1);
    p_i_3d = patientCutVers(idx_i,:);
    hole = Calc_Boundary(rootCutTris);
    idx_o = hole.boundary.edge(:,1);
    p_o_3d = rootCutVers(idx_o,:);
    proj = zeros(3,3);
    proj(:,3) = patientAxis';
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

%% 4.修整牙根的形状：

    % 牙冠和牙根融合
    finalVers = [patientCutVers; rootCutVers];
    reidx = [idx_i; idx_o + size(patientCutVers,1)];
    finalTris = [patientCutTris; rootCutTris + size(patientCutVers,1); reidx(Fa)];

    % 融合部分光滑
    A = adjacency_matrix(finalTris);
    L = A - diag(sparse(sum(A,2)));
    lhs = L;
    L2 = L*L;
    lhs(idx_o + size(patientCutVers,1), :) = L2(idx_o + size(patientCutVers,1), :);
    rhs = L*finalVers;
    rhs(idx_o + size(patientCutVers,1), :) = 0;
    finalVers = solve_equation(lhs, rhs, 1:size(patientCutVers,1), patientCutVers);
 
    
    writeOBJ('最终结果.obj', finalVers, finalTris);

    disp('finished');

