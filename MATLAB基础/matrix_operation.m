%% ��������

%% �����任

%% �߼�����

%% 
% ����
%       n = norm(X,p) ���ؾ��� X �� p-���������� p Ϊ 1��2 �� Inf��
%       ��� p = 1���� n �Ǿ������������֮�͡�
%       ��� p = 2���� n ������ max(svd(X))������ norm(X) ��Ч��
%       ��� p = Inf���� n �Ǿ������������֮�͡�
v1 = round(rand(1, 3) * 10);
m1 = round(rand(2, 3) * 10);
disp(m1);
disp(v1);
norm1 = norm(m1, 1);
disp(norm1);
disp(norm(v1, 1));
norm2 = norm(m1, 2);
disp(norm2);
disp(norm(v1, 2));