function handle = gui_repair_root(handle)
    load('D:\workstation\gitRepositories\matlabCode\tooth_root_test\matlab.mat')
    
    toothdata = handle.model.toothdata;
    nT = size(toothdata,2);
    
	% 牙齿中心和牙合中心
    toothCenter = zeros(nT,3);
    
    for i = 1:nT
        toothCenter(i,:) = mean(toothdata{1,i});
    end
    dentalCenter = mean(toothCenter);
    
    % 逐个牙齿进行拼接
    for i = 1:nT
        if i >= 5 && i <=10
            [Vt, Ft] = readOBJ('root1.obj');
        else
            [Vt, Ft] = readOBJ('root2.obj');
        end
        V0 = toothdata{1,i};
        F0 = toothdata{2,i};
        dir = toothdata{7,i};
        orig = toothdata{6,i};
        
        % 牙齿局部坐标系
        toothCoord(:,1) = normalizerow(cross(toothCenter(i,:) - dentalCenter, -dir))';
        toothCoord(:,3) = -dir';
        toothCoord(:,2) = cross(toothCoord(:,3)', toothCoord(:,1));
        
        % 提取待拼接牙齿的牙冠区域
        Vc = V0(1:toothdata{8,i},:);
        Fc = F0(sum(F0 <= toothdata{8,i}, 2)  == 3, :);
        
        % 提取边缘控制点
        n = 20;
        hole = Calc_Boundary(Fc);
        line = hole.boundary.edge(:,1);
        temp = ceil(linspace(1, length(line), n+1));
        pick = temp(1:end-1);
        cp0 = Vc(line(pick),:);
        
        % 径向基函数控制点
        eta = 0.1;
        cp = [bsxfun(@plus, cp0, eta*dir); cp0; bsxfun(@minus, cp0, eta*dir)];
        cv = [ones(n,1); zeros(n,1); -ones(n,1)];
        
        % 径向基函数
        dist2 = pdist2(cp, cp);
        P = [ones(3*n,1), cp];
        A = [dist2, P; P', zeros(4,4)];
        B = [cv; zeros(4,1)];
        coeff = A\B;
        
        % 标准牙xy缩放
        Vctmp = bsxfun(@minus, Vc, toothCenter(i,:)) * toothCoord;
        Vt(:,1) = Vt(:,1) * (max(Vctmp(:,1)) - min(Vctmp(:,1)))/(max(Vt(:,1)) - min(Vt(:,1)));
        Vt(:,2) = Vt(:,2) * (max(Vctmp(:,2)) - min(Vctmp(:,2)))/(max(Vt(:,2)) - min(Vt(:,2)));
        
        % 将标准牙的坐标系对齐
        Vt = bsxfun(@plus, toothCoord * Vt', toothCenter(i,:)')';
        
        % 画分割面
        syms x y z;
        f = coeff(1:3*n)' * ((x - cp(:,1)).^2 + (y - cp(:,2)).^2 + (z - cp(:,3)).^2).^(1/2) + ...
            [1 x y z] * coeff((3*n+1):end);
        f = matlabFunction(f);
        
        % 提取标准牙的牙根区域
        isRoot = f(Vt(:,1),Vt(:,2),Vt(:,3)) > 0;
        Fp = tt(Ft);
        tag = double( sum(isRoot(Ft),2) == 3 );
        C = connected_region(Fp, tag);
        [~, ind] = sort( cellfun('length', C), 'descend' );
        [Vr, Fr] = Remove_Point(Vt, Ft(C{ind(1)},:));
        
 
        % 画牙冠和牙根
        figure
        drawMesh(Vc, Fc, 'facecolor','y', 'edgecolor','none', 'facealpha',0.8);
        drawMesh(Vr, Fr, 'facecolor','g', 'edgecolor','none', 'facealpha',0.8);
        
        view(3)
        axis equal
        axis off
        camlight
        lighting gouraud
        set(gca, 'Position',[0 0 1 1]);
        
        % 在平面上建立连接关系
        hole = Calc_Boundary(Fc);
        idx_i = hole.boundary.edge(:,1);
        p_i_3d = Vc(idx_i,:);
        hole = Calc_Boundary(Fr);
        idx_o = hole.boundary.edge(:,1);
        p_o_3d = Vr(idx_o,:);
        proj = zeros(3,3);
        proj(:,3) = dir';
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
        
 
        % 牙冠和牙根融合
        V = [Vc; Vr];
        reidx = [idx_i; idx_o + size(Vc,1)];
        F = [Fc; Fr + size(Vc,1); reidx(Fa)];
        
        % 融合部分光滑
        A = adjacency_matrix(F);
        L = A - diag(sparse(sum(A,2)));
        lhs = L;
        L2 = L*L;
        lhs(idx_o + size(Vc,1), :) = L2(idx_o + size(Vc,1), :);
        rhs = L*V;
        rhs(idx_o + size(Vc,1), :) = 0;
        V = solve_equation(lhs, rhs, 1:size(Vc,1), Vc);
        
        figure
        drawMesh(V, F, 'facecolor','y', 'edgecolor','none', 'facealpha',1.0);
        
        view(3)
        axis equal
        axis off
        camlight
        lighting gouraud
        set(gca, 'Position',[0 0 1 1]);
  
        progressbar(i,nT);
    end
end