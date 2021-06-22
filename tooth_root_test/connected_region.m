function [C] = connected_region(TT, TAG)
    TOVISIT = 1;
    VISITED = 2;

    TAG(TAG ~= 0) = TOVISIT; % mark as to visit
    C = cell(1,1);
    Crow = 1;

    I = find(TAG);

    for i = I'
        if (TAG(i) == 1) % must be visited
            CElem = [i];
            toVisit = [i];
            TAG(i) = VISITED;

            while size(toVisit) ~= 0
                current = toVisit(1);
                toVisit = toVisit(2:end);
                
                neigh = TT(current,:);
                neigh = neigh(neigh ~= -1);
                neigh = neigh(TAG(neigh) == TOVISIT);
                
                TAG(neigh) = VISITED;
                CElem = [CElem, neigh];
                toVisit = [toVisit, neigh];
            end

            C{Crow} = CElem;
            Crow = Crow + 1;
        end
    end
end
