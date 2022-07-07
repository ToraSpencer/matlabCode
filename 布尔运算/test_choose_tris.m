%% ���������е����ཻ����Ƭ ���� ���д�Ľӿ�solve_self_intersection();
 clc;
 clear all;
%  [vers, tris] = readOBJ('./myData/output.obj');
%      outerTriIdx = 42505;                % �����ʷֺ������ǽ��沿�ֵ�һ������Ƭ��������ѡһ�����ɣ�
      [vers, tris] = readOBJ('./myData/output_from_cork.obj');
     outerTriIdx = 16992;                % �����ʷֺ������ǽ��沿�ֵ�һ������Ƭ��������ѡһ�����ɣ�
     versCount = size(vers, 1);
    trisCount = size(tris,1);

    
    % 1. ���������ڽӹ�ϵ
    edges = [tris(:,2) tris(:,3); tris(:,3) tris(:,1); tris(:,1) tris(:,2)];    % ����ߣ�������Ƭ�������Ӧ��
    edgesIdx = (1:3 * trisCount)';           % ����ߵ���������������
    
    Ei = sparse(edges(:,1), edges(:,2), 1);                    
    adjSM = Ei > 0;                                 % ���������ڽӾ����±�(i, j)��Ԫ��Ϊtrue == ���������(i, j);
    nonDlSM = adjSM - adjSM';            % ��˫�����Ϣ��  
    adjSM2 = adjSM + nonDlSM;         % 0��ʾ�ޱߣ�1��ʾ��˫��ߣ�2��-1��ʾ����ߣ�
    adjSM3 = Ei + Ei';                             % 0��ʾ�ޱߣ�1��ʾ����ߣ�2��ʾ˫��ߣ�
    
    % (i,j) > 0����������ij���ڣ�(i,j)��ֵ����������ij������������Ƭ����
     % (i,j) > 0����������ijΪ�߽�ߣ�(i,j)��ֵӦ��ֻ��1��-1
    [ii, jj, k1] = find(adjSM2);
    [ii, jj, k2] = find(adjSM3);        
 
    Ef = sparse(edges(:,1), edges(:,2), edgesIdx);          % ���������ڽӾ���Ȩ��Ϊ�ߵ�������
    SM1 = sparse(ii, jj, k1==1, versCount, versCount);
    SM2 = sparse(ii, jj, k2<=2, versCount, versCount);
    adj = ( SM1 & SM2).*Ef;    
    
    % �Ǳ߽粢�����αߣ�(i,j)��ֵ����������ij�ı��
     [ii, jj, si] = find(adj);                   % si��ֵ����ߵı��
    [ii, jj, v] = find(adj');                   % v��ֵ�������ߵı��
    
    E2F = repmat(1:trisCount,1,3)';      % �ɱߵı������������Ƭ������
    
    Fp = -ones(3*trisCount,1);      % ��ʼ�ڽӹ�ϵ����Ϊ-1
    Fp(si) = E2F(v);
    Fp = reshape(Fp, trisCount, 3); % ����Ƭÿ�����αߵ��ڽ�����Ƭ������
    
    [ii, jj] = find(Ei == 2);           % ���ҷ����α߼���
    [sel, k3] = ismember(edges, [ii, jj], 'rows');
    k3 = k3(sel);                       % ���ҷ����α��ڷ����α߼����е�λ��
    neid = edgesIdx(sel);                % ���ҷ����αߵı��
    [sel, k4] = ismember(edges, [jj, ii], 'rows'); % ���ҷ����α��ڽӵ�����Ƭ����
    k4 = k4(sel);                           % ���ҷ������ڽӱ��ڷ����α߼����е�λ��
    nfid = E2F(sel);                         % ���ҷ����αߵ��ڽ�����Ƭ������
    [useless, reidx] = sort(k4);
    nfid = reshape(nfid(reidx)', 2, [])';
    
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