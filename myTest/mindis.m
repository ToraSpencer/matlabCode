function[minValues,rows]=mindis(vers,compareVer, k)
% k == 1ʱ�� minValuesΪ����vers�еĵ㵽compareVer��������룬rowsΪ����������С�

versCount = length(vers);
temp = repmat(compareVer, versCount, 1);
temp = temp - vers;
temp = temp.^2;
Distance=sqrt(sum(temp,2));

x = [(1:length(Distance))' Distance];

[u,v]=sort(x(:,2));

minValues = u(1:k);

rows = v(1:k);

