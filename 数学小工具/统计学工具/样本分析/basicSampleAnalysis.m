%����һ��������������������ĸ���ͳ������
clc;
clear all;
X = [65 61 76 77 76 65 74 99 51 60 60 69 53 63 50 87 69 78 75 56 101 68];                       %����
n = numel(X);               %����������
X_bar = mean(X);            %������ֵ��
S = sqrt(1/(n-1)*(sum(X.*X)-n*X_bar^2));    %������׼�
chi = sum(X.*X);

%���������ݹ鷨��


%����λ����
X = sort(X);
if rem(n,2)==0
    median = (X(n/2)+X(n/2+1))/2;
else
    median = X(ceil(n/2));
end

mean = mean(X);

