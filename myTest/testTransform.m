clc
close all;
clear all;

%% 6.1 计算变形所需的参数——layers_from_handle()
mergedMesh = Read_Obj('合并网格.obj');
exterior = ReadObj('非合并区域的顶点');
mergedToothVers_copy = mergedMesh.vertex;         % 合并网格顶点数据
newTris = mergedMesh.face;                        % 合并网格三角片数据
vertex_count = size(mergedToothVers_copy, 1);       % 合并网格点数

Omega = 1:vertex_count;
temp = ismember(Omega,exterior);        % exterior三角片中的索引
temp = ~temp;                           
Omega = Omega(temp);        % 不在exterior三角片中的顶点索引  
 
shrinking_handle = exterior;
growing_interior = Omega;

interior_faces = intersect(newTris, newTris.*ismember(newTris,growing_interior), 'rows');
handle_faces = intersect(newTris, newTris.*~ismember(newTris,growing_interior), 'rows');
unionTris = union(handle_faces,interior_faces, 'rows');
H0 = setdiff(newTris, unionTris,'rows');
N0 = intersect(shrinking_handle,reshape(H0,1,size(H0,1)*size(H0,2)));