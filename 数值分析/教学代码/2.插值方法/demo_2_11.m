% ������ֵ����
clear; clc

X = [27.7, 28, 29, 30];
Y = [4.1, 4.3, 4.1, 3.0];
df0 = 3.0; dfn = -4.0;  % ��һ��߽�����

n = length(X) - 1;
H = diff(X); % ÿ����ֵС����ĳ���
mu = H(1:n-1) ./ (H(1:n-1) + H(2:n));
lambda = 1 - mu;
% ������ײ���
Y1 = diff(Y) ./ diff(X);  % һ�ײ���
Y2 = diff(Y1) ./ ( X(3:end) - X(1:end-2) );  % ���ײ���
% �����Ҷ���
d = 6*Y2;
d0 = 6/H(1) * (Y1(1) - df0);
dn = 6/H(end) * (dfn - Y1(end));
% ����ϵ������
A = 2*eye(n+1) + diag([mu(:);1],-1) +diag([1;lambda(:)],1);
b = [d0; d(:); dn];
M = A\b;

% ���� s_k(x) ��ϵ��, ������ʽ
p = zeros(n,4); p0 = zeros(n,4); 
for k = 1 : n
    p(k,1) = (M(k+1) - M(k))/(6*H(k));
    p(k,2) = M(k)/2;
    p(k,3) = (Y(k+1) - Y(k))/ H(k) - H(k)/6 * (2*M(k) + M(k+1));
    p(k,4) = Y(k);
    p0(k,1) = M(k)/(6*H(k));
    p0(k,2) = M(k+1)/(6*H(k));
    p0(k,3) = (Y(k) - M(k)*H(k)*H(k)/6) / H(k);
    p0(k,4) = (Y(k+1) - M(k+1)*H(k)*H(k)/6) / H(k);
end

% ���� Matlab ������������ֵ����
pp = spline(X, [df0;Y(:);dfn]);

%������
fprintf('����ʽһ���: \n'); disp(p0);
fprintf('����ʽ�����: \n'); disp(p);
fprintf('spline�Ľ��(��ʽ��): \n'); disp(pp.coefs);
