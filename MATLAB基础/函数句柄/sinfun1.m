function [y] = sinfun1(M)
x = 0:M-1;                      %没有预先开辟内存。
for k = 1:numel(x)
   y(k) = sin(x(k))/(100*pi); 
end
end

