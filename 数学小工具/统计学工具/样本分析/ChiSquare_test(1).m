%%
%正太总体标准差sigma的假设检验：
 clc;
 clear all;
 alpha = 0.01;
 sigma0 = 6.0;
%%
%(a)基本统计量直接给出：
n = 40;
S = 4.9;
%%
%(b)基本统计量需要计算：
X = [774 649 1210 546 431 612];
n = numel(X);
X_bar = mean(X);
S = sqrt(var(X));

%%
chi2 = (n-1)*S^2/(sigma0^2);
disp(strcat('检测统计量chi2 =',num2str(chi2)));
%%
%(a)算P值：左边检验：
p = chi2cdf(chi2,n-1);
disp(strcat('P值 =',num2str(p)));
%%
%(b)算P值：右边检验：
p = chi2cdf(t,n-1);
P = 1-p;
disp(strcat('P值 =',num2str(P)));
%%
%(c)算P值：双边检验：
p = tcdf(-abs(chi2),n-1);
P = 2*p;
disp(strcat('P值 =',num2str(P)));

%%
%(a)单边检验求t临界值：
beta = 1-alpha;
chi2_alpha = chi2inv(beta,n-1);
disp(strcat('chi2临界值为自由度为n-1的chi2分布上alpha分位点：',num2str(chi2_alpha)));
%%
%(a)双边检验求t临界值：
beta = 1-alpha/2;
chi2_HalfAlpha = chi2inv(beta,n-1);
disp(strcat('chi2临界值为自由度为n-1的chi2分布上alpha/2分位点：',num2str(chi2_HalfAlpha)));








%%
%卡方拟合检验
clc;
clear all;
k = 6;              %事件数。
e = [1/6 1/6 1/6 1/6 1/6 1/6];         %期望频率向量    
O = [27 31 42 40 28 32];                      %观测频数向量
alpha = 0.05;     %显著性水平

n = sum(O);          %样本容量
E = e*n;            %期望频数向量

chiSquare = sum((O-E).^2./E);
disp(strcat('检测统计量chiSquare =',num2str(chiSquare)));

beta = 1-alpha;
chiSquare_alpha = chi2inv(beta,k-1);
disp(strcat('卡方分布的上alpha分位点为',num2str(chiSquare_alpha)));

p = chi2cdf(chiSquare,k-1);
P = 1-p;
disp(strcat('P值 =',num2str(P)));








%%
%卡方独立性检验：
clc;
clear all;
r = 2;              
c = 2;
O_matrix = [14 1056;95 437];         %观测频数矩阵。
alpha = 0.05;     %显著性水平

O = reshape(O_matrix,1,r*c);      %观测频数向量,元素为按行取得。
n = sum(O);          %样本容量
rowSum = sum(O_matrix,2);
colSum = sum(O_matrix,1);
E_matrix = (rowSum*colSum)./n;
E = reshape(E_matrix,1,r*c);      %期望频率矩阵。

chiSquare = sum((O-E).^2./E);
disp(strcat('检测统计量chiSquare =',num2str(chiSquare)));

beta = 1-alpha;
chiSquare_alpha = chi2inv(beta,(r-1)*(c-1));
disp(strcat('卡方分布的上alpha分位点为',num2str(chiSquare_alpha)));

p = chi2cdf(chiSquare,(r-1)*(c-1));
P = 1-p;
disp(strcat('P值 =',num2str(P)));


%%
%卡方同质性检验：检测多个总体的同一分类的比例p是否相同。
clc;
clear all;
r = 6;             %总体数，总体在行上列出。
c = 1;              %变量可取值个数，个数在列上列出。
O_matrix = [27;31;42;40;28;32];         %观测频数矩阵。
alpha = 0.01;        %显著性水平


O = reshape(O_matrix,1,r*c);      %观测频数向量,元素为按行取得。
n = sum(O);          %样本容量
rowSum = sum(O_matrix,2);
colSum = sum(O_matrix,1);
E_matrix = (rowSum*colSum)./n;
E = reshape(E_matrix,1,r*c);

chiSquare = sum((O-E).^2./E);
disp(strcat('检测统计量chiSquare =',num2str(chiSquare)));

beta = 1-alpha;
chiSquare_alpha = chi2inv(beta,(r-1)*(c-1));
disp(strcat('卡方分布的上alpha分位点为',num2str(chiSquare_alpha)));

p = chi2cdf(chiSquare,(r-1)*(c-1));
P = 1-p;
disp(strcat('P值 =',num2str(P)));

