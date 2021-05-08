%�����������ľ�ֵ�Ƿ���ͬ��
clc;
clear all;
%%
%����������������ݡ�
k = 4;            %������
X1 = [1.21 0.57 0.56 0.13 1.3];   %����������
X2 = [0.94 0.87 0.46 0.58 1.03];
X3 = [0.07 0.66 0.1 0.82 0.94];
X4 = [0.85 1.78 1.47 2.25 1.64];
alpha = 0.05;     %������ˮƽ

X = [X1,X2,X3,X4];

n1 = numel(X1);
n2 = numel(X2);
n3 = numel(X3);
n4 = numel(X4);
n = n1+n2+n3+n4;

X1_bar = sum(X1)/n1;
X2_bar = sum(X2)/n2;
X3_bar = sum(X3)/n3;
X4_bar = sum(X4)/n4;
X_bar = sum(X)/n;

S1 = sqrt((sum(X1.^2)-n1*X1_bar^2)/(n1-1)); 
S2 = sqrt((sum(X2.^2)-n2*X2_bar^2)/(n2-1));
S3 = sqrt((sum(X3.^2)-n3*X3_bar^2)/(n3-1));
S4 = sqrt((sum(X4.^2)-n4*X4_bar^2)/(n4-1));
S = sqrt((sum(X.^2)-n*X_bar^2)/(n-1));
disp(strcat('����������ֵΪX_bar =  ',num2str(X_bar)));
disp(strcat('����������׼��ΪS = ',num2str(S)));
%%
%����SSTR������ƽ���͡�
SSTR = n1*(X1_bar-X_bar)^2+n2*(X2_bar-X_bar)^2+n3*(X3_bar-X_bar)^2+n4*(X4_bar-X_bar)^2;
disp(strcat('SSTR = ',num2str(SSTR)));

%%
%����SSE�����ƽ���ͣ�error sum of squares��
SSE = (n1-1)*S1^2+(n2-1)*S2^2+(n3-1)*S3^2+(n4-1)*S4^2;
disp(strcat('SSE = ',num2str(SSE)));

%%
%����SST/TSS����ƽ����(total sum of squares)
SST = SSTR + SSE;
disp(strcat('SST = ',num2str(SST)));

%%
%����MSTR���������
MSTR = SSTR/(k-1);
disp(strcat('MSTR = ',num2str(MSTR)));

%%
%����MSE:��������
MSE = SSE/(n-k);
disp(strcat('MSE = ',num2str(MSE)));

%%
%����F��������ֵ�Ĳ�����
F = MSTR/MSE;
disp(strcat('F = ',num2str(F)));

%%
%����Pֵ��
p = fcdf(F,k-1,n-k);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));

%%
%����F�ֲ�����alpha��λ�㣺
beta = 1-alpha;
F_alpha = finv(beta,k-1,n-k);
disp(strcat('���ɶ�Ϊ(k-1,n-k)��F�ֲ���alpha��λ��Ϊ��',num2str(F_alpha)));