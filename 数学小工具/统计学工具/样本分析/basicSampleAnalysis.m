%输入一组样本，输出该组样本的各种统计量。
clc;
clear all;
X = [65 61 76 77 76 65 74 99 51 60 60 69 53 63 50 87 69 78 75 56 101 68];                       %样本
n = numel(X);               %样本容量。
X_bar = mean(X);            %样本均值。
S = sqrt(1/(n-1)*(sum(X.*X)-n*X_bar^2));    %样本标准差。
chi = sum(X.*X);

%求众数：递归法：


%求中位数：
X = sort(X);
if rem(n,2)==0
    median = (X(n/2)+X(n/2+1))/2;
else
    median = X(ceil(n/2));
end

mean = mean(X);

