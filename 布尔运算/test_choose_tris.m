%% 查找网格中的自相交三角片 —— 斌杰写的接口solve_self_intersection();
 clc;
 clear all;
%  [vers, tris] = readOBJ('./myData/output.obj');
%      outerTriIdx = 42505;                % 三角剖分后的网格非交叉部分的一个三角片索引；任选一个即可；
      [vers, tris] = readOBJ('./myData/output_from_cork.obj');
     outerTriIdx = 16992;                % 三角剖分后的网格非交叉部分的一个三角片索引；任选一个即可；
     versCount = size(vers, 1);
    trisCount = size(tris,1);

    
    % 1. 构建网格邻接关系
    edges = [tris(:,2) tris(:,3); tris(:,3) tris(:,1); tris(:,1) tris(:,2)];    % 有向边，与三角片正法向对应；
    edgesIdx = (1:3 * trisCount)';           % 有向边的索引，列向量；
    
    Ei = sparse(edges(:,1), edges(:,2), 1);                    
    adjSM = Ei > 0;                                 % 顶点有向邻接矩阵；下标(i, j)的元素为true == 存在有向边(i, j);
    nonDlSM = adjSM - adjSM';            % 非双向边信息；  
    adjSM2 = adjSM + nonDlSM;         % 0表示无边，1表示是双向边，2和-1表示单向边；
    adjSM3 = Ei + Ei';                             % 0表示无边，1表示单向边，2表示双向边；
    
    % (i,j) > 0代表该有向边ij存在，(i,j)的值代表该有向边ij所关联的三角片数量
     % (i,j) > 0代表该有向边ij为边界边，(i,j)的值应该只有1或-1
    [ii, jj, k1] = find(adjSM2);
    [ii, jj, k2] = find(adjSM3);        
 
    Ef = sparse(edges(:,1), edges(:,2), edgesIdx);          % 顶点有向邻接矩阵，权重为边的索引；
    SM1 = sparse(ii, jj, k1==1, versCount, versCount);
    SM2 = sparse(ii, jj, k2<=2, versCount, versCount);
    adj = ( SM1 & SM2).*Ef;    
    
    % 非边界并且流形边，(i,j)的值代表该有向边ij的编号
     [ii, jj, si] = find(adj);                   % si的值代表边的编号
    [ii, jj, v] = find(adj');                   % v的值代表对面边的编号
    
    E2F = repmat(1:trisCount,1,3)';      % 由边的编号索引到三角片的索引
    
    Fp = -ones(3*trisCount,1);      % 初始邻接关系设置为-1
    Fp(si) = E2F(v);
    Fp = reshape(Fp, trisCount, 3); % 三角片每条流形边的邻接三角片的索引
    
    [ii, jj] = find(Ei == 2);           % 查找非流形边集合
    [sel, k3] = ismember(edges, [ii, jj], 'rows');
    k3 = k3(sel);                       % 查找非流形边在非流形边集合中的位置
    neid = edgesIdx(sel);                % 查找非流形边的编号
    [sel, k4] = ismember(edges, [jj, ii], 'rows'); % 查找非流形边邻接的三角片索引
    k4 = k4(sel);                           % 查找非流形邻接边在非流形边集合中的位置
    nfid = E2F(sel);                         % 查找非流形边的邻接三角片的索引
    [useless, reidx] = sort(k4);
    nfid = reshape(nfid(reidx)', 2, [])';
    
    Fn = -ones(trisCount,3,2);                 % 初始邻接关系设置为-1，每条非流形边对面存在2个三角片
    Fn([neid; neid+3*trisCount]) = [nfid(k3,1); nfid(k3,2)];

    % 2. 扩散查找连通区域
    TOVISIT = 1;        % 待访问
    VISITED = 2;        % 已访问
    
    TAG = repmat(TOVISIT, trisCount, 1); % 访问状态
    
    N = normalizerow(normals(vers, tris));
    
    toVisit = [outerTriIdx]; % 加入到待扩散序列中
    TAG(outerTriIdx) = VISITED;
    while ~isempty(toVisit)
        current = toVisit(1);
        toVisit = toVisit(2:end);
        for j = 1:3
            neigh = Fp(current,j);
            if (neigh == -1) % 判断当前三角片的边是否为流形边
                neigh = squeeze( Fn(current,j,:) );
                [useless, idx] = min( N(current,:) * N(neigh,:)' );
                neigh = neigh(idx);
            end

            if (TAG(neigh) == TOVISIT)
                TAG(neigh) = VISITED;
                toVisit = [toVisit neigh];
            end
        end
    end
    
    [versOut, trisOut] = Remove_Point(vers, tris(TAG==VISITED,:));
    
    writeOBJ('finalMesh.obj', versOut, trisOut);
    disp('finished.');