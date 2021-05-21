%% ��������
% �ҳ�/  ���\
clc;
clear all;
m1 = rand(3,3);
m2 = rand(3,3);
m3 = m1*m2;
temp1 = m3/m2;          % �ҳ��� m1 == m3/m2 == m3 * inv(m2);
temp2 = m1\m3;          % ����� m2 == m1\m3 == inv(m1)*m3;



%% �����任



%% �߼�����

%%    ����
%               n = norm(X,p) ���ؾ��� X �� p-���������� p Ϊ 1��2 �� Inf��
%               ��� p = 1���� n �Ǿ������������֮�͡�
%               ��� p = 2���� n ������ max(svd(X))������ norm(X) ��Ч��
%               ��� p = Inf���� n �Ǿ������������֮�͡�
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

%   �����󽻣�
A = [7 1 7 7 4]; 
B = [7 0 4 4 0];
[C,ia,ib] = intersect(A,B);
disp(C);            % ����AB�Ĺ���Ԫ��
disp(ia);           % ����Ԫ�ص�����
disp(ib);           % ����Ԫ�ص�����

A = [2 2 2; 0 0 1; 1 2 3; 1 1 1; 2, 2, 2];
B = [1 2 3; 2 2 2; 2 2 0];
[C,ia,ib] = intersect(A,B,'rows');
C1 = intersect(A,B,'rows');
disp(C);            % ����AB�Ĺ�����
disp(ia);           % �����е�����
disp(ib);           % �����е�����
disp(C1);


% ismember()���������������m1�е�ÿ��Ԫ����m2���Ƿ���ڡ�
m1 = [1,1,0,9,2,10, -1,1,15,1,3,1];
m2 = [1,2,3; 4, 5, 6; 7, 8, 9];
disp(ismember(m1, m2));

% union()��������Ĳ������Զ�ȥ���ظ�����
m1 = [1,2,3; 4,5,6; 1,2,3];
m2 = [-1,-2,-3; 4, 5, 6; 7, 8, 9; -4, -5, -6];
m3 = union(m1, m2, 'rows');
disp(m3);


% setdiff()��������Ԫ�صĲ�����Զ�ȥ���ظ�����
m1 = [1,2,3; 4,5,6; 1,2,3];
m2 = [-1,-2,-3; 4, 5, 6; 7, 8, 9; -4, -5, -6];
m3 = setdiff(m1, m2, 'rows');
disp(m3);



%% lu�ֽ⡪��[L,U,P,Q,D] = lu(S)�Ǹ�˹��Ԫ���ľ�����ʽ
% S == D * inv(P) * L * U * inv(Q)
% L�������Ǿ���U�������Ǿ���inv(P)������û�����inv(Q)���ҳ��û�����D�ǶԽ����ž���
clc;
clear all;
m1 = [1,2,3,4;5,6,7,8;9,10,0,12;13,14,15,0];
disp(rank(m1));
disp(det(m1));
sm1 = sparse(m1);
[sL, sU, sP, sQ, sD] = lu(sm1);
L = full(sL);
U = full(sU);
P = full(sP);
Q = full(sQ);
D = full(sD);
disp(m1);
disp(D * inv(P) * L * U * inv(Q) );
disp(D);



%% ϡ�����

%   sparse()�����ɳ��ܾ���õ�ϡ�����
%%
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

%%


%   full()������ϡ�����õ����ܾ���
m11 = full(sm1);
m2 = full(sm2);

% nnz()���� ����ϡ������еķ���Ԫ������
nonZero = nnz(sm1);
 
% nonzeros()���� ���� ϡ���������з���Ԫ�أ��洢��һ���������ڡ�
m3 = nonzeros(sm1);

% nzmax()���� ����Ϊϡ�����ķ��������Ĵ洢�ռ�����
size = nzmax(sm1);



%%
clc;
clear all;
m1 = [11,12,13,14;21,22,23,24; 31,32,33,34];
v1 = [1,2,3];
v2 = [3,4];
m2 = m1(v1, v2);
disp(m2);

V = v1'*v2;
disp(V);
