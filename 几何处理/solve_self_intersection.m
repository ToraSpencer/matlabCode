
% 斌杰写的网格自相交处理
function [versOut, trisOut] = solve_self_intersection(vers, tris, id)
    % 消除网格自相交
    trisCount = size(tris,1);
    
    % 1. 构建网格邻接关系
    E = [tris(:,2) tris(:,3); tris(:,3) tris(:,1); tris(:,1) tris(:,2)];
    
    Eid = (1:3 * trisCount)'; % 三角片边的编号
    Ef = sparse(E(:,1), E(:,2), Eid);
    Ei = sparse(E(:,1), E(:,2), 1);
    
    unadj = Ei > 0;                 % (i,j) > 0代表该有向边ij存在，(i,j)的值代表该有向边ij所关联的三角片数量
    nonadj = unadj - unadj'; % (i,j) > 0代表该有向边ij为边界边，(i,j)的值应该只有1或-1
    
    [ii, jj, k1] = find(unadj + nonadj);
    [ii, jj, k2] = find(Ei + Ei'); % (i,j)的值代表该无向边ij所关联的三角片数量
    adj = ( ...
        sparse(ii, jj, k1==1, size(Ef,1), size(Ef,2)) & ...
        sparse(ii, jj, k2<=2, size(Ef,1), size(Ef,2))).*Ef; % 非边界并且流形边，(i,j)的值代表该有向边ij的编号
    
    [ii, jj, si] = find(adj); % si的值代表边的编号
    [ii, jj, v] = find(adj'); % v的值代表对面边的编号
    
    E2F = repmat(1:trisCount,1,3)'; % 由边的编号索引到三角片的索引
    
    Fp = -ones(3*trisCount,1); % 初始邻接关系设置为-1
    Fp(si) = E2F(v);
    Fp = reshape(Fp, trisCount, 3); % 三角片每条流形边的邻接三角片的索引
    
    [ii, jj] = find(Ei == 2); % 查找非流形边集合
    [sel, k3] = ismember(E, [ii, jj], 'rows');
    k3 = k3(sel);               % 查找非流形边在非流形边集合中的位置
    neid = Eid(sel);            % 查找非流形边的编号
    [sel, k4] = ismember(E, [jj, ii], 'rows'); % 查找非流形边邻接的三角片索引
    k4 = k4(sel);           % 查找非流形邻接边在非流形边集合中的位置
    nfid = E2F(sel);            % 查找非流形边的邻接三角片的索引
    [useless, reidx] = sort(k4);
    nfid = reshape(nfid(reidx)', 2, [])';
    
    Fn = -ones(trisCount,3,2);                 % 初始邻接关系设置为-1，每条非流形边对面存在2个三角片
    Fn([neid; neid+3*trisCount]) = [nfid(k3,1); nfid(k3,2)];

    % 2. 扩散查找连通区域
    TOVISIT = 1;        % 待访问
    VISITED = 2;        % 已访问
    
    TAG = repmat(TOVISIT, trisCount, 1); % 访问状态
    
    N = normalizerow(normals(vers, tris));
    
    toVisit = [id]; % 加入到待扩散序列中
    TAG(id) = VISITED;
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
end
