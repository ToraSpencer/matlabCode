%% 查找网格中的自相交三角片 —— 斌杰写的接口solve_self_intersection();
     clc;
     clear all;
     
     debugFlag = 1;
    elemInfo.row = [];
    elemInfo.col = [];
    elemInfo.value = [];
     
    [vers, tris] = readOBJ('./myData/output_from_cork_test1.obj');
%     [vers, tris] = readOBJ('./myData/output_from_cork.obj');
    versCount = size(vers, 1);
    trisCount = size(tris,1);

    % 取非交叉部分的任意一个三角片——z坐标最大的顶点所在的一个三角片：
    zValues = vers(:, 3);
    [~, verIdx] = max(zValues);
    IM1 = (tris == verIdx);
    [rows, ~] = find(IM1);
    outerTriIdx = rows(1);
    
    if(debugFlag == 1)
        writeOBJ('inputMesh.obj', vers, tris);
        a = tris(outerTriIdx, 1);
        b = tris(outerTriIdx, 2);
        c = tris(outerTriIdx, 3);
        outerTriVers= [vers(a, :); vers(b, :); vers(c, :)]; 
        objWriteVertices('outerTriVers.obj', outerTriVers);
    end
    
    % 1. 构建网格邻接关系
    triEdges = [tris(:,2) tris(:,3); tris(:,3) tris(:,1); tris(:,1) tris(:,2)];    % 有向边，与三角片正法向对应；
    triEdgesIdx = (1:3 * trisCount)';                   % 有向边的索引，列向量；
    triEdgesCount = size(triEdgesIdx, 1);           % 有向边数；
    
    adjSM_double = sparse(triEdges(:,1), triEdges(:,2), 1); % 顶点有向邻接矩阵（double)；下标(i, j)的元素为n == 存在n条有向边(i, j);          
    adjSM = adjSM_double > 0;   % 顶点有向邻接矩阵（logical）；下标(i, j)的元素为1 == 存在有向边(i, j);
    
    % for debug
    [elemInfo.row, elemInfo.col, elemInfo.value] = find(adjSM_double);
    edgesCount = size(elemInfo.value, 1);      % 无向边个数（即两顶点之间若连接便视为存在一条无向边，至多一条）
    
    nonDlSM = adjSM - adjSM';            % 非双向边信息(double)；0表示无边或双向边；1表示ij边，-1表示ji边；
    adjSM2 = adjSM + nonDlSM;           % 0表示无边，1表示是双向边，2和-1表示单向边；
    adjSM3 = adjSM_double + adjSM_double';  % 0表示无边，(i, j)元素值表示ij边的数量（正反方向都算）；
    
    % for debug
    [elemInfo.row, elemInfo.col, elemInfo.value] = find(adjSM3);
 
    % (i,j) > 0代表该有向边ij存在，(i,j)的值代表该有向边ij所关联的三角片数量
     % (i,j) > 0代表该有向边ij为边界边，(i,j)的值应该只有1或-1
     left = [];
     right = [];
    [left, right, edgeLabel1] = find(adjSM2);    % 边方向列向量：0表示无边，1表示是双向边，2和-1表示单向边；
    [left, right, edgeLabel2] = find(adjSM3);    % 边重复数列向量：0表示无边，(i, j)元素值表示ij边的数量（正反方向都算）；
 
    wadjSM_double = sparse(triEdges(:,1), triEdges(:,2), triEdgesIdx);    % 顶点有向邻接矩阵，权重为边的索引；
    SM1 = sparse(left, right, edgeLabel1==1, versCount, versCount);     % 双向边的索引矩阵；
    SM2 = sparse(left, right, edgeLabel2<=2, versCount, versCount);    % 边重复数不大于2的索引矩阵；
    SM3 = SM1 & SM2;            % 无重复双向边（非边界且流形边）的索引矩阵；
    wMFadj_double = SM3.*wadjSM_double;    % 顶点流形边邻接矩阵(double)，权重为边的索引；
    
    % 非边界并且流形边，(i,j)的值代表该有向边ij的编号
     [left, right,  MFedgeIdx] = find(wMFadj_double);       % 流形边的索引；
    [left, right,  MFoppEdgeIdx] = find(wMFadj_double');    % 流形边对面边的索引（(a,b)边的对面边是(b, a)）
    
    % 三角片索引及其三条有向边索引的关系：索引为i的三角片三条有向边索引为：3*i-2, 3*i-1, 3*i;
    E2F = repmat(1:trisCount,1,3)';      % 由边的编号索引到三角片的索引    
    ttIdx = -ones(3*trisCount,1);                  % 三角片每条流形边的邻接三角片的索引，-1表示没有；
    ttIdx(MFedgeIdx) = E2F(MFoppEdgeIdx);
    ttIdx = reshape(ttIdx, trisCount, 3);         
    
    [left, right] = find(adjSM_double == 2);           % 查找非流形边集合
    [sel, k3] = ismember(triEdges, [left, right], 'rows');
    k3 = k3(sel);                                   % 查找非流形边在非流形边集合中的位置
    neid = triEdgesIdx(sel);                    % 查找非流形边的编号
    [sel, k4] = ismember(triEdges, [right, left], 'rows');  % 查找非流形边邻接的三角片索引
    k4 = k4(sel);                           % 查找非流形邻接边在非流形边集合中的位置
    nfid = E2F(sel);                         % 查找非流形边的邻接三角片的索引
    [useless, reidx] = sort(k4);
    temp = nfid(reidx)';
    nfid = reshape(temp, 2, [])';
    
    Fn = -ones(trisCount,3,2);                 % 初始邻接关系设置为-1，每条非流形边（交线段边）对面存在2个三角片
    Fn([neid; neid+3*trisCount]) = [nfid(k3,1); nfid(k3,2)];

    % 2. 扩散查找连通区域
    TOVISIT = 1;        % 待访问
    VISITED = 2;        % 已访问
    
    TAG = repmat(TOVISIT, trisCount, 1); % 访问状态
    
    N = normalizerow(normals(vers, tris));
    
    
    %% 三角片扩散：
    toVisit = [outerTriIdx];                        % 加入到待扩散序列中
    TAG(outerTriIdx) = VISITED;
    while ~isempty(toVisit)
        current = toVisit(1);
        toVisit = toVisit(2:end);
        for j = 1:3
            neigh = ttIdx(current,j);
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