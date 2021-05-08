%%
%��̫�����׼��sigma�ļ�����飺
 clc;
 clear all;
 alpha = 0.01;
 sigma0 = 6.0;
%%
%(a)����ͳ����ֱ�Ӹ�����
n = 40;
S = 4.9;
%%
%(b)����ͳ������Ҫ���㣺
X = [774 649 1210 546 431 612];
n = numel(X);
X_bar = mean(X);
S = sqrt(var(X));

%%
chi2 = (n-1)*S^2/(sigma0^2);
disp(strcat('���ͳ����chi2 =',num2str(chi2)));
%%
%(a)��Pֵ����߼��飺
p = chi2cdf(chi2,n-1);
disp(strcat('Pֵ =',num2str(p)));
%%
%(b)��Pֵ���ұ߼��飺
p = chi2cdf(t,n-1);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));
%%
%(c)��Pֵ��˫�߼��飺
p = tcdf(-abs(chi2),n-1);
P = 2*p;
disp(strcat('Pֵ =',num2str(P)));

%%
%(a)���߼�����t�ٽ�ֵ��
beta = 1-alpha;
chi2_alpha = chi2inv(beta,n-1);
disp(strcat('chi2�ٽ�ֵΪ���ɶ�Ϊn-1��chi2�ֲ���alpha��λ�㣺',num2str(chi2_alpha)));
%%
%(a)˫�߼�����t�ٽ�ֵ��
beta = 1-alpha/2;
chi2_HalfAlpha = chi2inv(beta,n-1);
disp(strcat('chi2�ٽ�ֵΪ���ɶ�Ϊn-1��chi2�ֲ���alpha/2��λ�㣺',num2str(chi2_HalfAlpha)));








%%
%������ϼ���
clc;
clear all;
k = 6;              %�¼�����
e = [1/6 1/6 1/6 1/6 1/6 1/6];         %����Ƶ������    
O = [27 31 42 40 28 32];                      %�۲�Ƶ������
alpha = 0.05;     %������ˮƽ

n = sum(O);          %��������
E = e*n;            %����Ƶ������

chiSquare = sum((O-E).^2./E);
disp(strcat('���ͳ����chiSquare =',num2str(chiSquare)));

beta = 1-alpha;
chiSquare_alpha = chi2inv(beta,k-1);
disp(strcat('�����ֲ�����alpha��λ��Ϊ',num2str(chiSquare_alpha)));

p = chi2cdf(chiSquare,k-1);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));








%%
%���������Լ��飺
clc;
clear all;
r = 2;              
c = 2;
O_matrix = [14 1056;95 437];         %�۲�Ƶ������
alpha = 0.05;     %������ˮƽ

O = reshape(O_matrix,1,r*c);      %�۲�Ƶ������,Ԫ��Ϊ����ȡ�á�
n = sum(O);          %��������
rowSum = sum(O_matrix,2);
colSum = sum(O_matrix,1);
E_matrix = (rowSum*colSum)./n;
E = reshape(E_matrix,1,r*c);      %����Ƶ�ʾ���

chiSquare = sum((O-E).^2./E);
disp(strcat('���ͳ����chiSquare =',num2str(chiSquare)));

beta = 1-alpha;
chiSquare_alpha = chi2inv(beta,(r-1)*(c-1));
disp(strcat('�����ֲ�����alpha��λ��Ϊ',num2str(chiSquare_alpha)));

p = chi2cdf(chiSquare,(r-1)*(c-1));
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));


%%
%����ͬ���Լ��飺����������ͬһ����ı���p�Ƿ���ͬ��
clc;
clear all;
r = 6;             %�������������������г���
c = 1;              %������ȡֵ�����������������г���
O_matrix = [27;31;42;40;28;32];         %�۲�Ƶ������
alpha = 0.01;        %������ˮƽ


O = reshape(O_matrix,1,r*c);      %�۲�Ƶ������,Ԫ��Ϊ����ȡ�á�
n = sum(O);          %��������
rowSum = sum(O_matrix,2);
colSum = sum(O_matrix,1);
E_matrix = (rowSum*colSum)./n;
E = reshape(E_matrix,1,r*c);

chiSquare = sum((O-E).^2./E);
disp(strcat('���ͳ����chiSquare =',num2str(chiSquare)));

beta = 1-alpha;
chiSquare_alpha = chi2inv(beta,(r-1)*(c-1));
disp(strcat('�����ֲ�����alpha��λ��Ϊ',num2str(chiSquare_alpha)));

p = chi2cdf(chiSquare,(r-1)*(c-1));
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));

