clc
close all;
clear all;

%% 6.1 �����������Ĳ�������layers_from_handle()
mergedMesh = Read_Obj('�ϲ�����.obj');
exterior = ReadObj('�Ǻϲ�����Ķ���');
mergedToothVers_copy = mergedMesh.vertex;         % �ϲ����񶥵�����
newTris = mergedMesh.face;                        % �ϲ���������Ƭ����
vertex_count = size(mergedToothVers_copy, 1);       % �ϲ��������

Omega = 1:vertex_count;
temp = ismember(Omega,exterior);        % exterior����Ƭ�е�����
temp = ~temp;                           
Omega = Omega(temp);        % ����exterior����Ƭ�еĶ�������  
 
shrinking_handle = exterior;
growing_interior = Omega;

interior_faces = intersect(newTris, newTris.*ismember(newTris,growing_interior), 'rows');
handle_faces = intersect(newTris, newTris.*~ismember(newTris,growing_interior), 'rows');
unionTris = union(handle_faces,interior_faces, 'rows');
H0 = setdiff(newTris, unionTris,'rows');
N0 = intersect(shrinking_handle,reshape(H0,1,size(H0,1)*size(H0,2)));