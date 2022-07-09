clc;
clear all;

%% 
%   sparse()�����ɳ��ܾ���õ�ϡ�����
clc;
clear all;
m1 = [ 0   0   0   5
      0   2   0   0
      1   3   0   0
      0   0   4   0];
sm1 = sparse(m1);
rowInfo = [1,2,3];
colInfo = [3,4,5];
valueList = [0.1, 0.2, 0.3];
sm2 = sparse(rowInfo, colInfo, valueList, 5, 6);

rowInfo = [1,2,3,4,5,6];
colInfo = [1,2,3,4,5,6];
valueList = [ 1 4
              2 5
              3 6 ];          % ȡԪ�ص�ʱ��������
sm3 = sparse(rowInfo, colInfo, valueList, 6, 6);

rowInfo = [1,1,3,4,5,5];
colInfo = [1,1,3,4,5,5];      % ��ͬ������Ϣ��ֵ�ۼ� 
valueList = [ 1 4
              2 5
              3 6 ];        
sm4 = sparse(rowInfo, colInfo, valueList, 6, 6);

%   full()������ϡ�����õ����ܾ���
m11 = full(sm1);
m2 = full(sm2);

% nnz()���� ����ϡ������еķ���Ԫ������
nonZero = nnz(sm1);
 
% nonzeros()���� ���� ϡ���������з���Ԫ�أ��洢��һ���������ڡ�
m3 = nonzeros(sm1);

% nzmax()���� ����Ϊϡ�����ķ��������Ĵ洢�ռ�����
size = nzmax(sm1);

% [row, col, value] = find(spM) �������ϡ������з���Ԫ�ص���Ϣ��
elemInfo.row = [];
elemInfo.col = [];
elemInfo.value = [];
[elemInfo.row, elemInfo.col, elemInfo.value] = find(sm1);

