function NodeVector = ChordLength_Para(Q,PreNodeVector,k)
N = size(Q,1);
n = N - 1;
NodeVector = zeros(1,n+k+2);
NodeVector(1, n+2 : n+k+2) = 1;
U = 0;
for i = k+2 :n+1
    for j = i-k:i-1
       U = U + PreNodeVector(1,j);
    end
    NodeVector(1,i) = U/k;
    U = 0;
end
