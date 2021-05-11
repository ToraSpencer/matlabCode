function[minValue,row]=matchest(C,P)
[N,M]=size(C);
Distance=zeros([1,N]);
Distance=sqrt(P.^2*ones(size(C'))+ones(size(P))*(C').^2-2*P*C');
[minValue,row]=min(Distance);