%  ƽ������ �ֽ�Ľ��շ�ʽ
function [L,x] = my_cholesky(A,b)

[m,n]=size(A);
if (m~=n)
    fprintf('\n Error: A ���Ƿ���!\n'); return; 
end

L = zeros(n); 

for k = 1:n
    L(k,k) = sqrt(A(k,k)-sum(L(k,1:k-1).^2));
    if (L(k,k) == 0)
        fprintf('\n Error: �� %d ����Ԫ��0 !\n',k); return; 
    end
    for j = k+1:n
        L(j,k) = (A(j,k) - sum( L(j,1:k-1).*L(k,1:k-1) ))/L(k,k);
    end
end

% ���Ly=b �� L'x=y
for i=1:n
    summary = 0;
    for t=1:i-1
        summary = summary + L(i,t) * y(t);
    end
    y(i)=(b(i) - summary) / L(i,i);
end
for i = n : -1 : 1
    summary = 0;
    for t=i+1:n
        summary = summary + L(t, i) * x(t);
    end
    x(i)=(y(i) - summary)/L(i,i);
end


% END
