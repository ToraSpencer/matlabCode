%  ����Ԫ Doolittle �ֽ�Ľ��շ�ʽ
function [L, U, x] = my_lu_with_column_pivoting(A, b)

[m,n]=size(A);
if (m~=n)
    fprintf('\n Error: A ���Ƿ���!\n'); return; 
end

L = eye(n); 
U = zeros(n);

Ip=1:n;  %���ĵ�k��Ԫ��M(k)��¼��k����Ԫ�����ڵ��к�

for k = 1:n
    % ���ֽ�PA=LU��P�����о��󣬼�ʵ��ѡ��Ԫ�أ��������н���
    [v, h] = max(abs(A(k:n,k)));
    h = k + h -1;
    Ip(k) = h ;
    temp = A(k,:);
    A(k,:) = A(h,:);
    A(h,:) = temp;
    
    A(k,k) = A(k,k)-sum(A(k,1:k-1)'.*A(1:k-1,k));
    if (A(k,k) == 0)
        fprintf('\n Error: �� %d ����Ԫ��0 !\n',k); return; 
    end
    for j = k+1:n
        % L��U����ֱ�Ӵ洢��A�����У����Խ�ʡ�洢�ռ�
        A(k,j) = A(k,j) - sum(A(k,1:k-1)'.*A(1:k-1,j));
        A(j,k) = (A(j,k) - sum( A(j,1:k-1)'.*A(1:k-1,k) ))/A(k,k);            
    end

end
L = tril(A,-1) + eye(n);
U = triu(A);

for i = 1:n
    t = Ip(i);
    if t ~= i
        temp = b(i);
        b(i) = b(t);
        b(t) = temp;
    end
end

y = 1:n;
y(1) = b(1);
for i=2:n
    summary = 0;
    for t=1:i-1
        summary = summary + L(i,t) * y(t);
    end
    y(i)=b(i) - summary;
    % ���շ�ʽ
%     t = (1 : i-1)';
%     y(i) = b(i) - sum(L(i,t) .* y(t));
end
x(n) = y(n)/U(n,n);
for i = n-1 : -1 : 1
    summary = 0;
    for t=i+1:n
        summary = summary + U(i,t) * x(t);
    end
    x(i)=(y(i) - summary)/U(i,i);
end

%END