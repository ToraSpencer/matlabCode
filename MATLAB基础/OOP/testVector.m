clc;
clear all;
%%
%�������ʵ����������һ������
vector0 = classSample_vector2D();
disp('����default constructor�����Ķ�ά��������');
disp(vector0);

vector1 = classSample_vector2D(1,2);
disp('����������ù��췽�������Ķ�ά��������');
disp(vector1);


vector1.normalize();
disp('���ù�һ������');
disp(vector1);

vector2 = classSample_vector3D(1,2,3);
disp('����������ù��췽����������ά��������');
disp(vector2);