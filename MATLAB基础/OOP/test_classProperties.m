clc;
clear all;
obj = N_randomNumber();

%%
a = obj.cumulativeProbability;      %����û����Ĭ��ֵ�����췽������û�и�ֵ����ʵ��������ĸ�������double�վ���
disp(class(a));

%%
%����һ��ʵ��������������Ծ������ˣ��̶����䣬����ʹ��dot���������set������Ϊ�ı䡣
%��������ÿ�α����ʶ�Ҫ���¼��㣬���������������к������������ʹ��ͬһ������ÿ�η��ʽ�����ܿ��ܲ�ͬ��
A = [];
B = [];
for i = 1:10
    A = [A;obj.value1];
    B = [B;obj.value2];
end
disp([A,B]);

%%
fprintf('��̬�ֲ��ı�׼��Ϊ��%d\n',N_randomNumber.sigma);   %Constant���Կ���ͨ����ֱ�ӵ��á�
fprintf('��̬�ֲ��ķ���Ϊ:%d\n',N_randomNumber.variance);   %һ��Constant���Կ��Դ�������һ��Constant���ԡ�

