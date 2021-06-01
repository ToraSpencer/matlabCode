clc;
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
 
%%  1.�и������

%       1.1 ȡ��������
x = 11;             % ȡ11����������
for i = 1:28
    if  dentalwithtooth1(i).ID == x
        rootTooth = dentalwithtooth1(i);
    end
end
 
for j = 1:28
    if (upax(j).fid == x)
        axisStandard = upax(j);
    end
end
for k = 1:28
    if (dental_crown(k).fid == x)
        crown = dental_crown(k);
    end
end

s = fix(x/10);  %ȡ��
g = mod(x,10);%ȡ��
 
temp = Read_Obj('��������׼��.obj');      % ԭ��������Ƭ���ˣ�������Լ������ġ�
rootTooth.vertices = temp.vertex;
rootTooth.faces = temp.face;

fdi = textread('FDIUpper__.dxt');
toothIdx = find (fdi == x);

namestr1 = ['toothUpper_',num2str(toothIdx-1),'.','obj'];
patientTooth = Read_Obj(namestr1);
namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];
gumline =  ReadObj(namestr2);
axis = ReadObj('AXISUpper_.obj');
axisPatient = axis(3*(toothIdx-1)+1:3*toothIdx,:);
crownTooth = crown.model;          % ��׼����


%      1.2 ȷ�����������и�㡣
centerPatient = mean(patientTooth.vertex);
y = max(gumline(:,2));
p_patient =[centerPatient(1),y(1),centerPatient(3)];
patientYdir = axisPatient(2,:);
cutPatientIdx = find(patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-2);
cutPatientCrownVers = patientTooth.vertex(cutPatientIdx,:);
    
    
%       1.3 ȷ�������и����ڱ߽硢����Ƭ��
a1 = ismember(patientTooth.face(:,1),cutPatientIdx);
aa1 = find(a1==1);
a2 = ismember(patientTooth.face(:,2),cutPatientIdx);
aa2 = find(a2==1);
a3 = ismember(patientTooth.face(:,3),cutPatientIdx);
aa3 = find(a3==1);
[c1,~,~] = intersect(aa1,aa2);
triIdxInCutPatient = intersect(c1,aa3); % ��������patientTooth��������������������cutPatientIdx���ҵ�������Ƭ��������

save('patientTooth.mat', 'patientTooth');
save('triIdxInCutPatient.mat', 'triIdxInCutPatient');


%�ұ߽�
trisInCutPatient = patientTooth.face(triIdxInCutPatient,:);
[raw_edges_list] = query_edges_list(trisInCutPatient,'sorted');
[~,~,iu] = unique(sort(raw_edges_list,2),'rows');
[i3,i4] = histc(iu,unique(iu));
lone_edges_idx_vect = i3(i4) == 1;
bdryEdges = unique(raw_edges_list(lone_edges_idx_vect,:),'rows');

%�߽�㣺
edgeVerIdx = unique([bdryEdges(:,1); bdryEdges(:,2)]);      % �������������߽������е�����������ظ�

%�߽�������ڵ��е�λ��
[~,tempRowIdx,~] = intersect(cutPatientIdx, edgeVerIdx);

%           ����������ԭ���������еĵ������� �����������и�֮��Ĳ��������еĵ�������
bdryEdges_newRep = zeros(size(bdryEdges));       % ����������ʾ�ı߽�ߡ�
for i = 1:length(tempRowIdx)
    index = find(bdryEdges == edgeVerIdx(i));
    bdryEdges_newRep(index) = tempRowIdx(i);
end

patientEdge = sotr_edge(bdryEdges_newRep,1);  % �����Ĳ����и����ڱ߽��,�����������ϵ����ų�(a, b);(b, c);(c, d);(d, e)��������ʽ
patientCutTris = ones(size(trisInCutPatient)); % ��������ʾ�Ĳ����и����������Ƭ��

for j = 1:length(cutPatientCrownVers)       
    index = find( trisInCutPatient == cutPatientIdx(j));
    patientCutTris(sub2ind(size(patientCutTris), index)) = j;
end

     
writeOBJ('�и��Ĳ�����������.obj', cutPatientCrownVers, patientCutTris)
%%
% 2. ���롪�����ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж���

%  2.1 ȷ����һ����תƽ���������
    centerCrown = mean(crownTooth.vertex);
    rootYdir = axisStandard.y;
    index_crown_biaozhun = find(crownTooth.vertex(:,2)<=(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(:,1) ...
                - centerCrown(1))+rootYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/rootYdir(2)));
    v_crown_biaozhun = abs(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(index_crown_biaozhun,1)-centerCrown(1))...
                 +rootYdir(3)*(crownTooth.vertex(index_crown_biaozhun,3)-centerCrown(3)))/rootYdir(2));
    in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
    point_crown =  crownTooth.vertex(index_crown_biaozhun(in_crown),:);

    p_tooth_biaozhun = mean(rootTooth.vertices);
    index_tooth_biaozhun = find(rootTooth.vertices(:,2)<=(p_tooth_biaozhun(2) - (rootYdir(1)*(rootTooth.vertices(:,1) ...
                     - p_tooth_biaozhun(1))+rootYdir(3)*(rootTooth.vertices(:,3) - p_tooth_biaozhun(3)))/rootYdir(2)));
    vz_tooth = abs(p_tooth_biaozhun(2) - (rootYdir(1)*(rootTooth.vertices(index_tooth_biaozhun,1)-p_tooth_biaozhun(1))...
               +rootYdir(3)*(rootTooth.vertices(index_tooth_biaozhun,3)-p_tooth_biaozhun(3)))/rootYdir(2));
    in_tooth = find(vz_tooth == max(vz_tooth));
    point_tooth = rootTooth.vertices(index_tooth_biaozhun(in_tooth),:);
    centcrownintooth = centerCrown+point_tooth-point_crown;     %�к��׼���ı�Ե����

    
%       2.2 ��һ����תƽ��
temp = inv([axisStandard.x;axisStandard.y;axisStandard.z]);
R = temp * axisPatient;
newCent = centcrownintooth*R;

moveVec = centerPatient - newCent;
movedRootTooth1 = (rootTooth.vertices) *R + repmat(moveVec,length(rootTooth.vertices),1);

save('movedRootTooth1.mat','movedRootTooth1');

%       2.3 ���icp�㷨��Ҫ�������㼯��
index_1 = find(patientTooth.vertex(:,2)>=(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))...
        & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+0.5);
point1 = patientTooth.vertex(index_1,:);
index_2 = find(movedRootTooth1(:,2)>=(centerPatient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2))...
        & movedRootTooth1(:,2)<=(p_patient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2))+0.5);
point2 = movedRootTooth1(index_2,:);



save('point1.mat','point1');
save('point2.mat','point2');

%       2.4 icp����
for i = 1:length(point2)
   [minValue,r]=mindis(point1,point2(i,:),1);
   row(i) = r;
end

point1 = point1(row,:);

[R,t,~,~,~,~] = icp(point1,point2);


%       2.5 �ڶ�����תƽ��
movedRootTooth2 = bsxfun(@plus,movedRootTooth1*R,t);


save('movedRootTooth2.mat', 'movedRootTooth2');
OBJwriteVertices('�ڶ�����תƽ�ƺ�Ĵ�����׼��.obj', movedRootTooth2);



%% 3. �и������׼��  

%       3.1 ȷ���и�㡣
rootCutIdx =  find(movedRootTooth2(:,2)>(p_patient(2) - ( patientYdir(1)*(movedRootTooth2(:,1) ...
                 - p_patient(1))+ patientYdir(3)*(movedRootTooth2(:,3) - p_patient(3)))/ patientYdir(2))-2);%�г���������������ֵ%ע�⣡���ղ������ڷ�����
rootCutVers = movedRootTooth2(rootCutIdx ,:);

a_t1 = ismember(rootTooth.faces(:,1),rootCutIdx);
aa_t1 = find(a_t1==1);
a_t2 = ismember(rootTooth.faces(:,2),rootCutIdx);
aa_t2 = find(a_t2==1);
a_t3 = ismember(rootTooth.faces(:,3),rootCutIdx);
aa_t3 = find(a_t3==1);
[c_t1,~,~] = intersect(aa_t1,aa_t2);
triIdxInCutRoot = intersect(c_t1,aa_t3);     % ��������root��������������������cutRoot���ҵ�������Ƭ��������

%       3.2 ȷ���и�����ı߽������Ƭ��
[raw_edges_list_t] = query_edges_list(rootTooth.faces(triIdxInCutRoot,:),'sorted');
[~,~,iu_t] = unique(sort(raw_edges_list_t,2),'rows');
[i3,i4] = histc(iu_t,unique(iu_t));
lone_edges_idx_vect_t = i3(i4) == 1;
rootBdryEdges = unique(raw_edges_list_t(lone_edges_idx_vect_t,:),'rows');       %  ����������ʾ���и������߽�ߡ�

trisInCutRoot = rootTooth.faces(triIdxInCutRoot,:);


%       3.3 ȷ���и���������Ƭ���߽�ӳ�䵽�ϲ������е�����Ƭ���߽硣
rootCutTrisInMergeMesh = ones (size(trisInCutRoot));           % ��������ʾ��cutRoot�е�����Ƭ��
for j = 1:length(rootCutVers)
    index = find( trisInCutRoot == rootCutIdx(j));
    rootCutTrisInMergeMesh(sub2ind(size(rootCutTrisInMergeMesh), index)) = j+length(cutPatientCrownVers);
end

%�߽�㣺
rootEdgeVersIdx = unique([rootBdryEdges(:,1);rootBdryEdges(:,2)]);  % �������������߽������е�����������ظ�

%�߽�����������е�λ��
[~,tempRowIdx,~] = intersect(rootCutIdx,rootEdgeVersIdx);

rootEdgeVers = rootCutVers(tempRowIdx,:);
tempRowIdx = tempRowIdx + length(cutPatientCrownVers);

rootBdryEdges_newRep = zeros(size(rootBdryEdges));          % ����������ʾ�ı߽�ߡ�
for i = 1:length(tempRowIdx)
    index = find(rootBdryEdges == rootEdgeVersIdx(i));
    rootBdryEdges_newRep(index) = tempRowIdx(i);
end


allTris = [patientCutTris; rootCutTrisInMergeMesh];            % �ϲ�����Ķ��㡣
allVers = [cutPatientCrownVers; rootCutVers];       % �ϲ���������Ƭ��δ�����


save('rootBdryEdges_newRep.mat', 'rootBdryEdges_newRep');
save('patientEdge.mat', 'patientEdge');     
save('allVers.mat', 'allVers');
save('allTris.mat', 'allTris');
save('rootEdgeVers.mat', 'rootEdgeVers');
save('bdryEdges_newRep.mat', 'bdryEdges_newRep');
writeOBJ('������Ƭǰ�ĺϲ�����.obj', allVers, allTris);

%% 4. �ϲ���������Ƭ

%5.20����
%1.���еı߽�㰴��һ����������
[~,r] = mindis(allVers(rootBdryEdges_newRep(:,1),:), allVers(bdryEdges_newRep(1,1),:), 1); 

% �ж�[edg_b(1,��),Edge_t(r,1)]de ������[Edge_t(r,��),edg_b(1,1)]�ķ����Ƿ���ͬ
% Ϊ�˷�ֹ���������������Edge_t(r,1)����ĵڶ�����
[r2,r3] = ind2sub(size(rootBdryEdges_newRep),find(rootBdryEdges_newRep == rootBdryEdges_newRep(r,2)));
rr = r2(ismember(r2,r)==0);
rm = r3(ismember(r2,r)==0);
tri11 = allVers(patientEdge(3,1),:) - allVers(patientEdge(1,1),:);
tri12 = allVers(rootBdryEdges_newRep(r,1),:) - allVers(patientEdge(3,1),:);
tri21 = allVers(rootBdryEdges_newRep(rr,3-rm),:) - allVers(rootBdryEdges_newRep(r,1),:);
tri22 = allVers(patientEdge(1,1),:) - allVers(rootBdryEdges_newRep(rr,3-rm),:);
s_b = sign(dot(allVers(patientEdge(1,1),:) -mean(rootEdgeVers),cross(tri11,tri12)));
s_t = sign(dot(allVers(rootBdryEdges_newRep(r,1),:) - mean(rootEdgeVers),cross(tri21,tri22)));

if s_b ~= s_t   % hit
    edg_t = sotr_edge(rootBdryEdges_newRep, r);
else
    E_t = [rootBdryEdges_newRep(:,2),rootBdryEdges_newRep(:,1)];
    edg_t = sotr_edge(E_t,r);
end

addTris = [];
for j = 1:length(edg_t)
    p1 = allVers(edg_t(j,1),:);
    p2 = allVers(edg_t(j,2),:);

    [~,ro]=mindis(allVers(patientEdge(:,1),:),(p1+p2)/2,1);
    roww(j) = ro;
    addTris = [addTris;[edg_t(j,:),patientEdge(ro,1)]];
end

if roww(1)<length(patientEdge)/2    % hit

    for j = 1:length(roww)-1
        if roww(j) ~= roww(j+1) 
            addTris = [addTris;[patientEdge(roww(j):roww(j+1)-1,:),repmat(edg_t(j,2),roww(j+1)-roww(j),1)]];
        end
    end

    if roww(end)==length(patientEdge)
        addTris = [addTris;[patientEdge(roww(end),:),edg_t(1,1)]];

        if roww(1)~=1
            addTris = [addTris;[patientEdge(1:roww(1),:),repmat(edg_t(1,1),roww(1),1)]];
        end

    elseif roww(end) > length(patientEdge)/2 && roww(end)<length(patientEdge)
        addTris = [addTris;[patientEdge(roww(end):length(patientEdge),:),repmat(edg_t(1,1),length(patientEdge)-roww(end),1)]];

        if roww(1)~=1
            addTris = [addTris;[patientEdge(1:roww(1),:),repmat(edg_t(1,1),roww(1),1)]];
        end

    elseif  roww(end) < length(patientEdge)/2 && roww(1)~=1

        s = find(roww > length(patientEdge)/2);
       addTris = [addTris;[patientEdge(roww(s(end)):length(patientEdge),:),repmat(edg_t(s(end),2),length(patientEdge)-roww(s(end))+1,1)]];
       addTris = [addTris;[patientEdge(1:roww(end),:),repmat(edg_t(s(end),2),roww(end),1)]]; 
       addTris = [addTris;[patientEdge(roww(end):roww(1)-1,:),repmat(edg_t(1,1),roww(1)-roww(end),1)]]; 
    end

else
   s = find(roww<length(patientEdge)/2);
   for j = s(1):length(roww)-1
        if roww(j) ~= roww(j+1) 
            addTris = [addTris;[patientEdge(roww(j):roww(j+1)-1,:),repmat(edg_t(j,2),roww(j+1)-roww(j),1)]];   
        end

   end 

   if s(1)>2
        for j = 1:s(1)-2
             if roww(j) ~= roww(j+1) 
                addTris = [addTris;[patientEdge(roww(j):roww(j+1)-1,:),repmat(edg_t(j,2),roww(j+1)-roww(j),1)]];   
             end
        end
   end

   if roww(end)<roww(1) 
      addTris = [addTris;[patientEdge(roww(end):roww(1)-1,:),repmat(edg_t(end,2),roww(1)-roww(end),1)]];
   end

   if roww(s(1)-1)<=length(patientEdge)
         addTris = [addTris;[patientEdge(roww(s(1)-1):length(patientEdge),:),repmat(edg_t(s(1)-1,2),length(patientEdge)-roww(s(1)-1)+1,1)]];
   end

   addTris = [addTris;[patientEdge(1:roww(s(1))-1,:),repmat(edg_t(s(1),1),roww(s(1))-1,1)]];

end

newTris = [allTris; addTris];

writeOBJ('����ǰ�ĺϲ�����.obj', allVers, newTris);
OBJwriteVertices('centerPatient.obj', centerPatient);

%% 5  ���ο��Ʋ���
index_crown_change = find(patientTooth.vertex(:,2)>(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-1 ...
                        & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+2);
patientTransform = patientTooth.vertex(index_crown_change,:);

%���������ֱ��γ�����
mergedToothVers = allVers;
index_tooth_change =  find(mergedToothVers(:,2)>(centerPatient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-1 ...
    &mergedToothVers(:,2)<=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+1);
mergeRegionVers = mergedToothVers(index_tooth_change,:);

% �������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�movedRootTooth2����crownforchange��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
indices = 1:size(mergedToothVers,1);

% �ϲ������У�����Ҫ������εĶ������������һ����������
exterior = indices(mergedToothVers(:,2)<(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-0.6...
    |mergedToothVers(:,2)>=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+3);
 

%% 6. ����
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
