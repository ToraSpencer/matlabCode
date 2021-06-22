function X = solve_equation(A, B, Acon, Bcon)
    % solve equation: AX = B with constraints that X(Acon) = Bcon
    if ~exist('Acon', 'var')
        Acon = [];
    end
    if ~exist('Bcon', 'var')
        Bcon = [];
    end
    
    if size(Acon,2) > 1
        Acon = Acon';
    end
    if length(Acon) ~= size(Bcon,1)
        Bcon = Bcon';
    end
    
%     [n, m] = size(A);    
%     r = ~ismember((1:n)', Acon);
%     
%     lhs = [A(r,:);sparse(1:length(Acon), Acon, 1, length(Acon), m)];
%     rhs = [B(r,:);Bcon];
%     
%     X = lhs\rhs;

    [n, m] = size(A);

    lhs = A;
    rhs = B;
    
    lhs(Acon,:) = sparse(1:length(Acon), Acon, 1, length(Acon), m);
    rhs(Acon,:) = Bcon;
    
    X = lhs\rhs;
end