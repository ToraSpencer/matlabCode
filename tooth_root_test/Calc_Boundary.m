function hole = Calc_Boundary(face)
    % 	compute the boundary of a 3D mesh(边的顺序与原始边界三角片的边顺序相反)
    %
    %   hole = Calc_Boundary(face);

    E = [face(:,2) face(:,3); ...
         face(:,3) face(:,1); ...
         face(:,1) face(:,2)];
    S = sparse(E(:,1), E(:,2), 1);
    S = S - S';
    [i,j,s] = find(S == -1);
    be = [i,j];
    if isempty(be)
        hole = [];
        return;
    end
    
    nhole = 0;
    unseen = true(length(be), 1);
    while any(unseen)
        nhole = nhole + 1;
        
        idx = find(unseen, 1);
        start = be(idx, 1);
        hole.boundary(nhole).edge(1,:) = be(idx,:);
        unseen(idx) = false;        
        next = be(idx, 2);
        
        while next ~= start
            idx = (be(:,1) == next & unseen);
            hole.boundary(nhole).edge = [hole.boundary(nhole).edge; be(idx,:)];
            unseen(idx,:) = false;    
            next = be(idx, 2);
        end

        hole.boundary(nhole).nBoundary = size(hole.boundary(nhole).edge, 1);
    end

    hole.nHole = nhole;
end