clc;
clear all;
%!!!���ݱ�����������
alpha = 0.05;

X = [102 101 94 79 79]';
Y = [175 169 182 146 144]';
n = numel(X);
X1_bar = mean(X);
X2_bar = mean(Y);
S1 = sqrt(var(X));
S2 = sqrt(var(Y));

%%
%�����������ϵ��
r = corr(X,Y,'type','Pearson');
disp(strcat('�������ϵ��r =',num2str(r)));

%%
%��С���˷����
%�ع鷽��Ϊy = beta1*x+beta0; beta1��б�ʣ�beta0�ǽؾࡣ
%b1,b0�ֱ����������Ƶ�б�ʺͽؾࡣ
b1 = (n*sum(X.*Y)-sum(X)*sum(Y))/(n*sum(X.^2)-sum(X)^2);
b0 = mean(Y)-b1*mean(X);
disp(strcat('��Ϸ���б�ʹ���ֵΪb1 =',num2str(b1)));
disp(strcat('��Ϸ��̽ؾ�Ϊb0 =',num2str(b0)));

%%
%���������ϵ����=0�ļ�����飺
%����裺��ȫ���Բ���أ���=0
t = abs(r)/sqrt((1-r^2)/(n-2));      %���ɶ�Ϊn-2
disp(strcat('���ͳ����t =',num2str(t)));

%%
%����Pֵ
p = tcdf(t,n-2);
P = 1-p;
disp(strcat('Pֵ =',num2str(P)));

%%
%t���ٽ�ֵ��
beta = 1-alpha/2;
t_HalfAlpha = tinv(beta,n-2);
disp(strcat('t���ٽ�ֵΪt�ֲ���alpha/2��λ�㣺',num2str(t_HalfAlpha)));

%%
%����ˮƽalpha��r���ٽ�ֵ��
r_critical = abs(t_HalfAlpha)/sqrt(t_HalfAlpha^2+n-2)





%%
%�в������
beta1 = 1;
beta0 = 2;
e = Y-(beta1*X+beta0);
SSE = sum(e.^2);
n = numel(e);
Se = sqrt(SSE/(n-2));





%%
%��������ع�����б���Ƿ�Ϊ0������x��ֵ����Ԥ��yֵ�Ƿ����������б��Ϊ0��yֵ��һ������������Ҫx��Ԥ�⣩
%�����:б��Ϊ0��beta1 == 0 
clc;
clear all;
t = (b1-beta1)/(Se/sqrt(Sxx));      %�������ɶ�Ϊn-2��t�ֲ�



