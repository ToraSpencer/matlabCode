function[minValues,rows]=mindis(patientTransform,vers,k)

% [N,M]=size(C);
% Distance=zeros([1,N]);
Distance=sqrt(sum((repmat((vers),length(patientTransform),1) - patientTransform).^2,2));
x = [(1:length(Distance))' Distance];
[u,v]=sort(x(:,2));
minValues = u(1:k);
rows = v(1:k);

