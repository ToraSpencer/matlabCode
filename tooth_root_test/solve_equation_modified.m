function X = solve_equation_modified(A, B, Acon, Bcon)
    % solve equation: AX = B with constraints that X(Acon) = Bcon

    Acon = Acon';
 
    [~, m] = size(A);
    temp = sparse(1:length(Acon), Acon, 1, length(Acon), m);
    A(Acon,:) = temp;
    B(Acon,:) = Bcon;
    
    X = A\B;
end