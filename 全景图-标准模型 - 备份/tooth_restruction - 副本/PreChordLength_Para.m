function PreNodeVector = PreChordLength_Para(Q)
[N M]= size(Q);
n = N -1;
PreNodeVector = zeros(1, n+1);
PreNodeVector(1,n+1) = 1;
Len = zeros(1, n-1);
if(M == 2)
    for iP = 2: N
    Len(iP-1) = sqrt((Q(iP,1) - Q(iP-1,1))^2 ...
        + (Q(iP,2) - Q(iP-1,2))^2);
%         + (Q(iP,3) - Q(iP-1,3))^2);
    end
end
if(M == 3)
    for iP = 2: N
    Len(iP-1) = sqrt((Q(iP,1) - Q(iP-1,1))^2 ...
        + (Q(iP,2) - Q(iP-1,2))^2 ...
        + (Q(iP,3) - Q(iP-1,3))^2);
    end
end
Lsum = sum(Len);
for iP = 2: n
    PreNodeVector(1, iP) = PreNodeVector(1, iP-1) + Len(iP-1)/Lsum;
end



