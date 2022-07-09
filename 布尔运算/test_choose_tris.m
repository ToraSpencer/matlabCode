%% ���������е����ཻ����Ƭ ���� ���д�Ľӿ�solve_self_intersection();
     clc;
     clear all;
     
    elemInfo.row = [];
    elemInfo.col = [];
    elemInfo.value = [];
     
%     [vers, tris] = readOBJ('./myData/output_from_cork_test1.obj');
    [vers, tris] = readOBJ('./myData/output_from_cork.obj');
    versCount = size(vers, 1);
    trisCount = size(tris,1);
    
    % ȡ�ǽ��沿�ֵ�����һ������Ƭ����z�������Ķ������ڵ�һ������Ƭ��
    zValues = vers(:, 3);
    [~, verIdx] = max(zValues);
    IM1 = (tris == verIdx);
    [rows, ~] = find(IM1);
    outerTriIdx = rows(1);
    
    % 1. ���������ڽӹ�ϵ
    triEdges = [tris(:,2) tris(:,3); tris(:,3) tris(:,1); tris(:,1) tris(:,2)];    % ����ߣ�������Ƭ�������Ӧ��
    triEdgesIdx = (1:3 * trisCount)';           % ����ߵ���������������
    
    adjSM_double = sparse(triEdges(:,1), triEdges(:,2), 1); % ���������ڽӾ���double)���±�(i, j)��Ԫ��Ϊn == ����n�������(i, j);          
    adjSM = adjSM_double > 0;   % ���������ڽӾ���logical�����±�(i, j)��Ԫ��Ϊ1 == ���������(i, j);
    
    % for debug
    [elemInfo.row, elemInfo.col, elemInfo.value] = find(adjSM_double);
    edgesCount = size(elemInfo.value, 1);      % ����߸�����
    
    nonDlSM = adjSM - adjSM';            % ��˫�����Ϣ(double)��0��ʾ�ޱ߻�˫��ߣ�1��ʾij�ߣ�-1��ʾji�ߣ�
    adjSM2 = adjSM + nonDlSM;         % 0��ʾ�ޱߣ�1��ʾ��˫��ߣ�2��-1��ʾ����ߣ�
    adjSM3 = adjSM_double + adjSM_double';  % 0��ʾ�ޱߣ�(i, j)Ԫ��ֵ��ʾij�ߵ����������������㣩��
    
    % for debug
    [elemInfo.row, elemInfo.col, elemInfo.value] = find(adjSM3);
 
    % (i,j) > 0����������ij���ڣ�(i,j)��ֵ����������ij������������Ƭ����
     % (i,j) > 0����������ijΪ�߽�ߣ�(i,j)��ֵӦ��ֻ��1��-1
     left = [];
     right = [];
    [left, right, edgeLabel1] = find(adjSM2);    % �߷�����������0��ʾ�ޱߣ�1��ʾ��˫��ߣ�2��-1��ʾ����ߣ�
    [left, right, edgeLabel2] = find(adjSM3);    % ���ظ�����������0��ʾ�ޱߣ�(i, j)Ԫ��ֵ��ʾij�ߵ����������������㣩��
 
    wadjSM_double = sparse(triEdges(:,1), triEdges(:,2), triEdgesIdx);    % ���������ڽӾ���Ȩ��Ϊ�ߵ�������
    SM1 = sparse(left, right, edgeLabel1==1, versCount, versCount);     % ˫��ߵ���������
    SM2 = sparse(left, right, edgeLabel2<=2, versCount, versCount);    % ���ظ���������2����������
    SM3 = SM1 & SM2;            % ���ظ�˫��ߣ��Ǳ߽������αߣ�����������
    wMFadj_double = SM3.*wadjSM_double;    % �������α��ڽӾ���(double)��Ȩ��Ϊ�ߵ�������
    
    % �Ǳ߽粢�����αߣ�(i,j)��ֵ����������ij�ı��
     [left, right,  MFedgeIdx] = find(wMFadj_double);       % si��ֵ����ߵı��
    [left, right,  v] = find(wMFadj_double');                   % v��ֵ�������ߵı��
    
    E2F = repmat(1:trisCount,1,3)';      % �ɱߵı������������Ƭ������
    
    Fp = -ones(3*trisCount,1);      % ��ʼ�ڽӹ�ϵ����Ϊ-1
    Fp(MFedgeIdx) = E2F(v);
    Fp = reshape(Fp, trisCount, 3); % ����Ƭÿ�����αߵ��ڽ�����Ƭ������
    
    [left, right] = find(adjSM_double == 2);           % ���ҷ����α߼���
    [sel, k3] = ismember(triEdges, [left, right], 'rows');
    k3 = k3(sel);                       % ���ҷ����α��ڷ����α߼����е�λ��
    neid = triEdgesIdx(sel);                % ���ҷ����αߵı��
    [sel, k4] = ismember(triEdges, [right, left], 'rows'); % ���ҷ����α��ڽӵ�����Ƭ����
    k4 = k4(sel);                           % ���ҷ������ڽӱ��ڷ����α߼����е�λ��
    nfid = E2F(sel);                         % ���ҷ����αߵ��ڽ�����Ƭ������
    [useless, reidx] = sort(k4);
    temp = nfid(reidx)';
    nfid = reshape(temp, 2, [])';
    
    Fn = -ones(trisCount,3,2);                 % ��ʼ�ڽӹ�ϵ����Ϊ-1��ÿ�������α߶������2������Ƭ
    Fn([neid; neid+3*trisCount]) = [nfid(k3,1); nfid(k3,2)];

    % 2. ��ɢ������ͨ����
    TOVISIT = 1;        % ������
    VISITED = 2;        % �ѷ���
    
    TAG = repmat(TOVISIT, trisCount, 1); % ����״̬
    
    N = normalizerow(normals(vers, tris));
    
    toVisit = [outerTriIdx]; % ���뵽����ɢ������
    TAG(outerTriIdx) = VISITED;
    while ~isempty(toVisit)
        current = toVisit(1);
        toVisit = toVisit(2:end);
        for j = 1:3
            neigh = Fp(current,j);
            if (neigh == -1) % �жϵ�ǰ����Ƭ�ı��Ƿ�Ϊ���α�
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