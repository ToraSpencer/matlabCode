function CoefMatrix = GetControlMatrix(Q,PreNodeVector,NodeVector)
[M N] = size(Q);
Niu = zeros(M-2,M);
for i = 1:M-2
    for j = 1:M
        Niu(i,j) = BaseFunction(j-1,3,PreNodeVector(i+1),NodeVector);
    end
end
coef_matrix = eye(M);
coef_matrix(2:M-1,:) = Niu;
inverse_matrix = pinv(coef_matrix);
CoefMatrix = inverse_matrix*Q;





