%% ��������
% �ҳ�/  ���\
clc;
clear all;
m1 = rand(3,3);
m2 = rand(3,3);
m3 = m1*m2;
temp1 = m3/m2;          % �ҳ��� m1 == m3/m2 == m3 * inv(m2);
temp2 = m1\m3;          % ����� m2 == m1\m3 == inv(m1)*m3;


%% �������ݵ���ɾ��ģ�

% min(), max()����������ֵ��
clc;
clear all;
vec1 = [8.1, 7.3, 6.2, 5.5, 4.4, 3.1, 2.5, 1.1];
[value, index] = min(vec1);
 
% [row, col, value] = find(M) ������þ����з���Ԫ�ص���Ϣ��
M1 =  [0, 8.1, 7.3; 6.2, 0, 4.4; -3.1, 2.5, 1.1];
elemInfo.row = [];
elemInfo.col = [];
elemInfo.value = [];
[elemInfo.row, elemInfo.col, elemInfo.value] = find(M1);

% [row, col, value] = find(Ԫ������) ���� ��þ���������������Ԫ�ص���Ϣ��
% [row, col] = find(Ԫ������)  
[elemInfo.row, elemInfo.col, elemInfo.value] = find(M1 < 0);        % �����value���߼�ֵ1


%% ʹ����������Է����飺
% 3x+y = 1; 2x-y = 2; 
clc;
clear all;
A = [3,1; 2,-1];
b = [1;2];              % A*x = b;
x1 = inv(A)*b;
x2 = A\b;               % A���������ϡ�����
disp(x1);
disp(x2);



%% �����任



%% �߼�����



%%    ����
%               n = norm(X), ���ؾ���X��2��������sqrt(x1^2+x2^2+x3^2 ....);
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

%   �����󽻣�intersect()����
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



%% ����ֵ�ֽ⡪��svd()
clc;
clear all;

%       [U,S,V] = svd(A) ִ�о��� A ������ֵ�ֽ⣬��� A = U*S*V'������S�ǶԽǾ���Ԫ��ΪA������ֵ��
A = [1 0 1; -1 -2 0; 0 1 -1];
[U,S,V] = svd(A); 
disp(U);
disp(S);
disp(V);

%       s = svd(A) �Խ���˳�򷵻ؾ��� A ������ֵ��


%% ��������
clc;
clear all;

% ������������ʽ����һ������������������������Ԫ��Ϊ1���������Ϊ0��
M1 =  [0, 8.1, 7.3; 6.2, 0, 4.4; -3.1, 2.5, 1.1];
IM1 = (M1 < 0);     
M2 = M1.*IM1;               % ���߼������ˣ�����ɸѡԪ�أ�

 
