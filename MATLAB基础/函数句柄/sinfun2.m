function [y] = sinfun2(M)
x = 0:M-1;
y = zeros(1,numel(x));              %Ԥ�ȿ����ڴ�
for k = 1:numel(x)
   y(k) = sin(x(k))/(100*pi); 
end
end

