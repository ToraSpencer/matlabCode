% Newton ǰ�幫ʽ
clear;
h = 0.1;
X = 0.0 : h : 0.5;
% Y = [1.00000,0.99500,0.98007,0.95534,0.92106,0.87758];
Y = cos(X);

% ������ײ��
n = length(X) - 1;
p = zeros(n+1,n+1);
p(:,1) = Y(:);
for k = 1 : n
    p(1:n-k+1,k+1) = diff(p(1:n-k+2,k));
end

% ���� 4 �� Newton ǰ�幫ʽ
q = p(1,:);
x = 0.048; t = (x-X(1))/h;
y = q(1); z = 1;
for k = 1 : 4 % ������ 4������ n(=5)
    z = z*(t-k+1)/k;
    y = y + z*q(k+1);
end
fprintf('�Ĵ�Newtonǰ�幫ʽ�ļ�����Ϊ: cos(%.3f)=%.5f\n',x,y);
