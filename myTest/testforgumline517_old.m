clc
clear all
functionname='testforgumline517.m'; functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir 'trianglerayintersection'])
addpath([functiondir 'toothmodel'])
addpath([functiondir 'theratofp'])
addpath([functiondir 'Root deformation'])
addpath([functiondir 'modelread'])
addpath([functiondir 'mixFE'])
addpath([functiondir 'margin line'])
addpath([functiondir 'Arch curve'])
addpath([functiondir 'Dentalmodelwithaxis'])
addpath([functiondir '������루��׼���Ͳ������ڶ��룩'])
addpath([functiondir 'siofmodeltooth'])
addpath([functiondir 'attachments'])
addpath([functiondir 'modelStadard'])
addpath([functiondir 'testdata'])
addpath([functiondir 'mesh process'])
addpath([functiondir 'testdata(4)'])
% % ���ر�׼��
load('dental_crown.mat');
load('dentalmodelwithroot0.1forlow.mat');
load('axisofdental_crowm');
load('trfromtoothtocrown.mat');

%%  1.
x = 11;             % ȡ11����������
for i = 1:28
    if  dentalwithtooth1(i).ID == x
        rootTooth = dentalwithtooth1(i);         % �������ı�׼������
    end
end
 
for j = 1:28
    if (upax(j).fid == x)
        axisStandard = upax(j);
    end
end

for k = 1:28
    if (dental_crown(k).fid == x)
        crown = dental_crown(k);            % ֻ�����ڵı�׼�����񣿣���
    end
end


temp = Read_Obj('��������׼��.obj');      % ԭ��������Ƭ���ˣ�������Լ������ġ�
rootTooth.vertices = temp.vertex;
rootTooth.faces = temp.face;


s = fix(x/10);      %ȡ��
g = mod(x,10);      %ȡ��
 

fdi = textread('FDIUpper__.dxt');
toothIdx = find (fdi == x);
 
namestr1 = ['toothUpper_',num2str(toothIdx-1),'.','obj'];
patientTooth = Read_Obj(namestr1);
namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];
gumline =  ReadObj(namestr2);
axis = ReadObj('AXISUpper_.obj');
axisPatient = axis(3*(toothIdx-1)+1:3*toothIdx,:);
crownTooth = crown.model;          % ��׼����



%       1.2 ȷ���������ڵ��и��
centerPatient = mean(patientTooth.vertex);
y = max(gumline(:,2));
p_patient =[centerPatient(1),y(1),centerPatient(3)];
patientYdir = axisPatient(2,:);


cutPatientIdx = find(patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-2);
cutPatientCrownVers = patientTooth.vertex(cutPatientIdx,:);


%       1.3 ȷ���������ڱ߽�
a1 = ismember(patientTooth.face(:,1),cutPatientIdx);
aa1 = find(a1==1);      
a2 = ismember(patientTooth.face(:,2),cutPatientIdx);
aa2 = find(a2==1);
a3 = ismember(patientTooth.face(:,3),cutPatientIdx);
aa3 = find(a3==1);
[c1,~,~] = intersect(aa1,aa2);
triIdxInCutPatient = intersect(c1,aa3);     % ��������patientTooth��������������������cutPatientIdx���ҵ�������Ƭ��������


save('patientTooth.mat', 'patientTooth');
save('triIdxInCutPatient.mat', 'triIdxInCutPatient');


%�ұ߽�
trisInCutPatient = patientTooth.face(triIdxInCutPatient,:);
[raw_edges_list] = query_edges_list(trisInCutPatient,'sorted');


[~,~,iu] = unique(raw_edges_list,'rows');           % ����iuΪ�������ε������������ȱȽϵ�һ��Ԫ�أ��ٱȽϵڶ���Ԫ�ء�
[i3,i4] = histc(iu,unique(iu));
lone_edges_idx_vect = i3(i4) == 1;          

temp = raw_edges_list(lone_edges_idx_vect,:);       % i3 == i4����
bdryEdges = unique(temp,'rows');

save('raw_edges_list.mat', 'raw_edges_list');

%�߽�㣺
edgeVerIdx = unique([bdryEdges(:,1); bdryEdges(:,2)]);      % �������������߽������е�����������ظ�

%�߽�������ڵ��е�λ��
[~,tempRowIdx,~] = intersect(cutPatientIdx, edgeVerIdx);
edgeVers = cutPatientCrownVers(tempRowIdx,:);            % ���б߽��

%           ����������ԭ���������еĵ������� �����������и�֮��Ĳ��������еĵ�������

bdryEdges_newRep = zeros(size(bdryEdges));       % ����������ʾ�ı߽�ߡ�
for i = 1:length(tempRowIdx)
    index = find(bdryEdges == edgeVerIdx(i));
    bdryEdges_newRep(index) = tempRowIdx(i);
end


patientEdge = sotr_edge(bdryEdges_newRep,1);     % �����Ĳ����и����ڱ߽��,�����������ϵ����ų�(a, b);(b, c);(c, d);(d, e)��������ʽ
patientCutTris = ones(size(trisInCutPatient));  % ��������ʾ�Ĳ����и����������Ƭ��

for j = 1:length(cutPatientCrownVers)       
    index = find( trisInCutPatient == cutPatientIdx(j));
    patientCutTris(sub2ind(size(patientCutTris), index)) = j;
end


%        1.4 ȷ��������׼���и����
centerCrown = mean(crownTooth.vertex);
rootYdir = axisStandard.y;
index_crown_biaozhun = find(crownTooth.vertex(:,2)<=(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(:,1) ...
            - centerCrown(1))+rootYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/rootYdir(2)));
v_crown_biaozhun = abs(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(index_crown_biaozhun,1)-centerCrown(1))...
             +rootYdir(3)*(crownTooth.vertex(index_crown_biaozhun,3)-centerCrown(3)))/rootYdir(2));
in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
point_crown =  crownTooth.vertex(index_crown_biaozhun(in_crown),:);

%��׼����
p_tooth_biaozhun = mean(rootTooth.vertices);
index_tooth_biaozhun = find(rootTooth.vertices(:,2)<=(p_tooth_biaozhun(2) - (rootYdir(1)*(rootTooth.vertices(:,1) ...
                 - p_tooth_biaozhun(1))+rootYdir(3)*(rootTooth.vertices(:,3) - p_tooth_biaozhun(3)))/rootYdir(2)));
vz_tooth = abs(p_tooth_biaozhun(2) - (rootYdir(1)*(rootTooth.vertices(index_tooth_biaozhun,1)-p_tooth_biaozhun(1))...
           +rootYdir(3)*(rootTooth.vertices(index_tooth_biaozhun,3)-p_tooth_biaozhun(3)))/rootYdir(2));
in_tooth = find(vz_tooth == max(vz_tooth));
point_tooth = rootTooth.vertices(index_tooth_biaozhun(in_tooth),:);
centcrownintooth = centerCrown+point_tooth-point_crown;   %�к��׼���ı�Ե����


%%
% 2. ���롪�����ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж���
temp = inv([axisStandard.x;axisStandard.y;axisStandard.z]);
R = temp * axisPatient;
newCent = centcrownintooth*R;

%       2.1 ��һ����תƽ��
moveVec = centerPatient - newCent;
movedRootTooth1 = (rootTooth.vertices) *R + repmat(moveVec,length(rootTooth.vertices),1);

%����仯���ݵ���̬
index_1 = find(patientTooth.vertex(:,2)>=(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))...
        & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+0.5);
point1 = patientTooth.vertex(index_1,:);        % ���������������������£�y�����ϸ���p_patient�ĵ㼯
index_2 = find(movedRootTooth1(:,2)>=(centerPatient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2))...
        & movedRootTooth1(:,2)<=(p_patient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2))+0.5);
point2 = movedRootTooth1(index_2,:);            % ��һ����תƽ�ƺ�������ڲ�����������ϵ�£�y�����ϸ���p_patient�ĵ㼯
 

save('movedRootTooth1.mat', 'movedRootTooth1');
save('point1.mat', 'point1');
save('point2.mat', 'point2');
writeOBJ('��һ����תƽ�ƺ�Ĵ�����׼��.obj',movedRootTooth1, rootTooth.faces)
OBJwriteVertices('point1.obj', point1);
OBJwriteVertices('point2.obj', point2);

%       2.2 �ڶ�����תƽ��
for i = 1:length(point2)
   [minValue,r]=mindis(point1,point2(i,:),1);
   row(i) = r;
end

point1 = point1(row,:);

[R,t,BRt,e,~,~] = icp(point1,point2);
movedRootTooth2 = bsxfun(@plus,movedRootTooth1*R,t);

save('movedRootTooth2.mat', 'movedRootTooth2');
OBJwriteVertices('�ڶ�����תƽ�ƺ�Ĵ�����׼��.obj', movedRootTooth2);


%%
% 3. �и��ȷ��������׼�����и�֡�
rootCutIdx =  find(movedRootTooth2(:,2)>(p_patient(2) - ( patientYdir(1)*(movedRootTooth2(:,1) ...
                 - p_patient(1))+ patientYdir(3)*(movedRootTooth2(:,3) - p_patient(3)))/ patientYdir(2))-2);%�г���������������ֵ%ע�⣡���ղ������ڷ�����
rootCutVers = movedRootTooth2(rootCutIdx ,:);

OBJwriteVertices('�ϲ�������������ֵ㼯.obj', rootCutVers);



%  ��������Ƭ��Ϣ
a_t1 = ismember(rootTooth.faces(:,1),rootCutIdx);
aa_t1 = find(a_t1==1);
a_t2 = ismember(rootTooth.faces(:,2),rootCutIdx);
aa_t2 = find(a_t2==1);
a_t3 = ismember(rootTooth.faces(:,3),rootCutIdx);
aa_t3 = find(a_t3==1);
[c_t1,~,~] = intersect(aa_t1,aa_t2);
c_t2 = intersect(c_t1,aa_t3);   % ��������root��������������������cutRoot���ҵ�������Ƭ��������

%�ұ߽�
[raw_edges_list_t] = query_edges_list(rootTooth.faces(c_t2,:),'sorted');
[~,~,iu_t] = unique(sort(raw_edges_list_t,2),'rows');
[i3,i4] = histc(iu_t,unique(iu_t));
lone_edges_idx_vect_t = i3(i4) == 1;
lone_edges_list_t = unique(raw_edges_list_t(lone_edges_idx_vect_t,:),'rows');

Edge_t = zeros(size(lone_edges_list_t));

 
ff_t = rootTooth.faces(c_t2,:);
f_t0 = ones (size(ff_t));
for j = 1:length(rootCutVers)
    index = find( ff_t == rootCutIdx(j));
    f_t0(sub2ind(size(f_t0), index)) = j+length(cutPatientCrownVers);
end
f_end = [patientCutTris;f_t0];v_end = [cutPatientCrownVers;rootCutVers];

%�߽�㣺
edge_p_t_local_list = unique([lone_edges_list_t(:,1);lone_edges_list_t(:,2)]);

%�߽�����������е�λ��
[C_t,edge_p_t_local,~] = intersect(rootCutIdx,edge_p_t_local_list);

edge_p_t = rootCutVers(edge_p_t_local,:);
edge_p_t_local = edge_p_t_local+length(cutPatientCrownVers);
Edge_t = zeros(size(lone_edges_list_t));
for i = 1:length(edge_p_t_local)
    index = find(lone_edges_list_t == edge_p_t_local_list(i));
    Edge_t(index) = edge_p_t_local(i);
end


%5.20����
%1.���еı߽�㰴��һ����������
[~,r]=mindis(v_end(Edge_t(:,1),:),v_end(bdryEdges_newRep(1,1),:),1);

%�ж�[edg_b(1,��),Edge_t(r,1)]de ������[Edge_t(r,��),edg_b(1,1)]�ķ����Ƿ���ͬ
%Ϊ�˷�ֹ���������������Edge_t(r,1)����ĵڶ�����
[r2,r3] = ind2sub(size(Edge_t),find(Edge_t == Edge_t(r,2)));
rr = r2(ismember(r2,r)==0);rm = r3(ismember(r2,r)==0);
tri11 = v_end(patientEdge(3,1),:) - v_end(patientEdge(1,1),:);tri12 = v_end(Edge_t(r,1),:) - v_end(patientEdge(3,1),:);
tri21 = v_end(Edge_t(rr,3-rm),:) - v_end(Edge_t(r,1),:);tri22 = v_end(patientEdge(1,1),:) - v_end(Edge_t(rr,3-rm),:);
s_b = sign(dot(v_end(patientEdge(1,1),:) -mean(edgeVers),cross(tri11,tri12)));
s_t = sign(dot(v_end(Edge_t(r,1),:) - mean(edge_p_t),cross(tri21,tri22)));

if s_b ~= s_t
   edg_t = sotr_edge(Edge_t,r);
else
    E_t = [Edge_t(:,2),Edge_t(:,1)];
    edg_t = sotr_edge(E_t,r);
end


dt = [];

for j = 1:length(edg_t)
    p1 = v_end(edg_t(j,1),:);
    p2 = v_end(edg_t(j,2),:);

    [~,ro]=mindis(v_end(patientEdge(:,1),:),(p1+p2)/2,1);
    roww(j) = ro;

    dt = [dt;[edg_t(j,:),patientEdge(ro,1)]];


end
if roww(1)<length(patientEdge)/2
    for j = 1:length(roww)-1
        if roww(j) ~= roww(j+1) 
            dt = [dt;[patientEdge(roww(j):roww(j+1)-1,:),repmat(edg_t(j,2),roww(j+1)-roww(j),1)]];
        end
    end
    if roww(end) > length(patientEdge)/2 && roww(end)<length(patientEdge)
        dt = [dt;[patientEdge(roww(end):length(patientEdge),:),repmat(edg_t(1,1),length(patientEdge)-roww(end),1)]];
    elseif  roww(end) < length(patientEdge)/2
       dt = [dt;[patientEdge(1:roww(end),:),repmat(edg_t(end,1),roww(end),1)]]; 
    end
    if roww(1)~=1
        dt = [dt;[patientEdge(1:roww(1)-1,:),repmat(edg_t(1,1),roww(1)-1,1)]];
    end
else
   s = find(roww<length(patientEdge)/2);
   roww = [roww roww(1:s(1)-1)];
   for j = s(1):length(roww)-1
        if roww(j) ~= roww(j+1) 
            dt = [dt;[patientEdge(roww(j):roww(j+1),:),repmat(edg_t(j,2),roww(j+1)-roww(j)+1,1)]];   
        end

   end 
   dt = [dt;[patientEdge(1:roww(s(1))-1,:),repmat(edg_t(1,2),roww(s(1))-1,1)]];

end
b = select_holes_and_boundary(v_end,[f_end;dt]);

pats = [];



%% 5  ���ο��Ʋ���
index_crown_change = find(patientTooth.vertex(:,2)>(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-1 ...
                        & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+2);
patientTransform = patientTooth.vertex(index_crown_change,:);

%���������ֱ��γ�����
mergedToothVers = v_end;
index_tooth_change =  find(mergedToothVers(:,2)>(centerPatient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-1 ...
    &mergedToothVers(:,2)<=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+1);
mergeRegionVers = mergedToothVers(index_tooth_change,:);



%�������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�tooth_T����crownforchange��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
indices = 1:size(mergedToothVers,1);
%�����ε�����
exterior = indices(mergedToothVers(:,2)<(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-0.6...
    |mergedToothVers(:,2)>=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+3);
 
   
%% 6. ����
newTris = [f_end;dt];
[omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(mergedToothVers,1), newTris, exterior);
[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(mergedToothVers,newTris, 'ext','voronoi', 'no_flatten',omega,N0,N1);

%�ҵ��������ָ����ڱ��β�������ĵ㣬�������ϵĵ��������ڴ�
finalVers = mergedToothVers;
for i = 1:length(mergeRegionVers)
   [minValue,r]=mindis(patientTransform,mergeRegionVers(i,:),1);
   minvalue(i) = minValue;  row(i) = r;
   finalVers(index_tooth_change(i),:) = patientTransform(r,:);
end

BZ1 = zeros(size(mergedToothVers,1),3);
finalVers = biharm_solve_with_factor( ...
    bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
    newTris, finalVers, omega, N0, N1, 'ext',  'no_flatten',BZ1,mergedToothVers);
 
 
writeOBJ('��������.obj',finalVers, newTris)


disp('finished.');
