%  ��˹��ȥ��
function [x] = my_ge_with_column_pivoting(A, b)

[m,n]=size(A);
if (m~=n)
    fprintf('\n Error: A ���Ƿ���!\n'); return; 
end

m = zeros(n);

for k = 1:n-1
    % Ѱ������Ԫ
  [v, h] = max(abs(A(k:n,k)));
  h = k + h -1;
  %�����У�����b
  if k ~= h
    temp = A(k,:);
    A(k,:) = A(h,:);
    A(h,:) = temp;
    temp = b(k);
    b(k) = b(h);
    b(h) = temp;  
  end
  %��Ԫ����
  for i = k+1:n
    m(i,k) = A(i,k)/A(k,k);
    A(i,:) = A(i,:) - m(i,k) * A(k,:);
    b(i) = b(i) - m(i,k) * b(k);
  end
end
%�ش�����
x(n) = b(n)/A(n,n);
for k = n-1 : -1 : 1
    summary = 0;
    for j = k+1 : n
        summary = summary + A(k,j) * x(j);
    end
    x(k) = (b(k) - summary) / A(k,k);
end
        