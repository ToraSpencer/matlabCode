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
 
 
x = 41;

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
crownTooth = crown.model;          % ��׼����
        
        
s = fix(x/10);  %ȡ��
g = mod(x,10);%ȡ��

    
if s == 1 || s == 2
    fdi = textread('FDIUpper__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
        disp( '�ò���û�д����� ');
        return;
    else
        fdi = textread('FDIUpper__.dxt');
        toothIdx = find (fdi == x);
        namestr1 = ['toothUpper_',num2str(toothIdx-1),'.','obj'];
        patientTooth = Read_Obj(namestr1);
        namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);
        axis = ReadObj('AXISUpper_.obj');
        axisPatient = axis(3*(toothIdx-1)+1:3*toothIdx,:);

   
axisStandardVers = [axisStandard.x; axisStandard.y; axisStandard.z];
OBJwriteVertices('axisPatient.obj', axisPatient);
OBJwriteVertices('axisStandard.obj', axisStandardVers);
OBJwriteVertices('gumline.obj', gumline);
writeOBJ('patientTooth.obj',patientTooth.vertex, patientTooth.face);
writeOBJ('crownTooth.obj',crownTooth.vertex, crownTooth.face);
writeOBJ('rootTooth.obj',rootTooth.vertices, rootTooth.faces);

%%  1.�и������

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
        triIdxInCutPatient = intersect(c1,aa3);
        trisInCutPatient = patientTooth.face(triIdxInCutPatient,:);
        
        [raw_edges_list] = query_edges_list(patientTooth.face(triIdxInCutPatient,:),'sorted');
        [~,~,iu] = unique(sort(raw_edges_list,2),'rows');
        [i3,i4] = histc(iu,unique(iu));
        lone_edges_idx_vect = i3(i4) == 1;
        bdryEdges = unique(raw_edges_list(lone_edges_idx_vect,:),'rows');
        
        %�߽�㣺
        edgeVerIdx = unique([bdryEdges(:,1);bdryEdges(:,2)]);
        
        %�߽�������ڵ��е�λ��
        [~,tempRowIdx,~] = intersect(cutPatientIdx,edgeVerIdx);
        
 %           ����������ԭ���������еĵ������� �����������и�֮��Ĳ��������еĵ�������
        edge_p_b = cutPatientCrownVers(tempRowIdx,:);
        patientBdryEdges_newRep = zeros(size(bdryEdges));
        for i = 1:length(tempRowIdx)
            index = find(bdryEdges == edgeVerIdx(i));
            patientBdryEdges_newRep(index) = tempRowIdx(i);
        end
        
        patientEdge = sotr_edge(patientBdryEdges_newRep,1);    
        patientCutTris = ones (size(trisInCutPatient)); % ��������ʾ�Ĳ����и����������Ƭ��
        
        for j = 1:length(cutPatientCrownVers)
            index = find( trisInCutPatient == cutPatientIdx(j));
            patientCutTris(sub2ind(size(patientCutTris), index)) = j;
        end
        
        
writeOBJ('�и��Ĳ�����������.obj', cutPatientCrownVers, patientCutTris)



%% 2. ���롪�����ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж���

%  2.1 ȷ����һ����תƽ���������
        centerCrown = mean(crownTooth.vertex);
        rootYdir = axisStandard.y;
        tempCrownIdx = find(crownTooth.vertex(:,2)<=(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(:,1) ...
                    - centerCrown(1))+rootYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/rootYdir(2)));
        crownValue = abs(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(tempCrownIdx,1)-centerCrown(1))...
                     +rootYdir(3)*(crownTooth.vertex(tempCrownIdx,3)-centerCrown(3)))/rootYdir(2));
        maxIdx = find(crownValue == max(crownValue));
        point_crown =  crownTooth.vertex(tempCrownIdx(maxIdx),:);

        %��׼����
        centerRoot = mean(rootTooth.vertices);
        tempRootIdx = find(rootTooth.vertices(:,2)<=(centerRoot(2) - (rootYdir(1)*(rootTooth.vertices(:,1) ...
                         - centerRoot(1))+rootYdir(3)*(rootTooth.vertices(:,3) - centerRoot(3)))/rootYdir(2)));
        rootValue = abs(centerRoot(2) - (rootYdir(1)*(rootTooth.vertices(tempRootIdx,1)-centerRoot(1))...
                   +rootYdir(3)*(rootTooth.vertices(tempRootIdx,3)-centerRoot(3)))/rootYdir(2));
        maxIdx = find(rootValue == max(rootValue));
        point_tooth = rootTooth.vertices(tempRootIdx(maxIdx),:);
        centcrownintooth = centerCrown+point_tooth-point_crown;%�к��׼���ı�Ե����
        
        
        %%���ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж��룬������ղ������ڵ��и�λ�ý����и�
        temp = inv([axisStandard.x;axisStandard.y;axisStandard.z]);
        R = temp * axisPatient;
        movingVec = centerPatient - centcrownintooth*R;
        movingMat = repmat(movingVec,length(rootTooth.vertices),1);
        movedRootTooth1 = (rootTooth.vertices) *R + movingMat;
        
        OBJwriteVertices('��һ����תƽ�ƺ�Ĵ�����׼��.obj', movedRootTooth1);


        %����仯���ݵ���̬
        index_1 = find(patientTooth.vertex(:,2)>=(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))...
                & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+0.5);
        point1 = patientTooth.vertex(index_1,:);
        
        index_2 = find(movedRootTooth1(:,2)>=(centerPatient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2))...
                & movedRootTooth1(:,2)<=(p_patient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2))+0.5);
        point2 = movedRootTooth1(index_2,:);
        
        
%       2.4 icp����
        for i = 1:length(point2)
           [minValue,startEdgeIdx]=mindis(point1,point2(i,:),1); 
           row(i) = startEdgeIdx;
        end
        
        point11 = point1(row,:);
        [R,t,~,~,~,~] = icp(point11,point2);
        
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

writeOBJ('������Ƭǰ�ĺϲ�����.obj', allVers, allTris);
OBJwriteVertices('�и���������.obj', rootCutVers);
writeOBJ('�и���������.obj', rootCutVers, trisInCutRoot);


%% 4. �ϲ���������Ƭ

%5.20����
%1.���еı߽�㰴��һ����������
re1 = patientEdge(1,1);                             % �������ڱ�Ե���һ�����������
rp11 = allVers(re1,:);
rbeVersIdx = rootBdryEdges_newRep(:,1);                         % ������Ե�������
[~,startEdgeIdx] = mindis(allVers(rbeVersIdx,:), rp11, 1);       % ������Ե���о������ڱ�Ե��һ����������������

% �ж�[patientEdge(1,��),Edge_t(r,1)]de ������[Edge_t(r,��),patientEdge(1,1)]�ķ����Ƿ���ͬ
% Ϊ�˷�ֹ���������������Edge_t(r,1)����ĵڶ�����
afterStartIdx = rootBdryEdges_newRep(startEdgeIdx,2);       % ��ʼ�ߵĺ󶥵��������
[row1,col1] = ind2sub(size(rootBdryEdges_newRep),find(rootBdryEdges_newRep == afterStartIdx));     

rr = row1(ismember(row1,startEdgeIdx)==0);      % row1����Ԫ���в�����startIdx�е��Ǹ�
rm = col1(ismember(row1,startEdgeIdx)==0);      % 1����2������ѡȡ��������Ӧ��λ�á�

pe3 = patientEdge(3,1);
pe1 = patientEdge(1,1);
pp1 = allVers(pe1,:);
pp3 = allVers(pe3,:) ;
rp1 = allVers(rootBdryEdges_newRep(startEdgeIdx,1),:);
rp3 = allVers(rootBdryEdges_newRep(rr,3-rm),:);
rCenter = mean(rootEdgeVers);

tri11 = pp3 - pp1;
tri12 = rp1 - pp3;
tri21 = rp3 - rp1;
tri22 = pp1 - rp3;

s_b = dot(pp1 - rCenter, cross(tri11,tri12));
s_t = dot(rp1 - rCenter, cross(tri21,tri22));

if s_b * s_t <= 0   % 11hit  ��ͬ�� 
    rootEdges = sotr_edge(rootBdryEdges_newRep, startEdgeIdx);
else            % 12hit
    E_t = [rootBdryEdges_newRep(:,2),rootBdryEdges_newRep(:,1)];
    rootEdges = sotr_edge(E_t,startEdgeIdx);
end


%   4.2 
addTris = [];
middle = [];


for j = 1:length(rootEdges)
    pe1 = allVers(rootEdges(j,1),:);
    p2 = allVers(rootEdges(j,2),:);
    middle(j, :) = (pe1+p2)/2;
    [~,ro]=mindis(allVers(patientEdge(:,1),:),(pe1+p2)/2,1);
    roww(j) = ro;
    addTris = [addTris;[rootEdges(j,:),patientEdge(ro,1)]];
end

OBJwriteVertices('middle.obj',middle);


%   4.3. 
if roww(1) < length(patientEdge)/2    % 11hit

    for j = 1:length(roww)-1
        if roww(j) ~= roww(j+1) 
            addTris = [addTris;[patientEdge(roww(j):roww(j+1)-1,:),repmat(rootEdges(j,2),roww(j+1)-roww(j),1)]];
        end
    end

    if roww(end) == length(patientEdge)
        addTris = [addTris; [patientEdge(roww(end),:), rootEdges(1,1)]];

        if roww(1)~=1
            addTris = [addTris;[patientEdge(1:roww(1),:),repmat(rootEdges(1,1),roww(1),1)]];
        end

    elseif roww(end) > length(patientEdge)/2 && roww(end)<length(patientEdge)
        addTris = [addTris;[patientEdge(roww(end):length(patientEdge),:),repmat(rootEdges(1,1),length(patientEdge)-roww(end)+1,1)]];        
        if roww(1)~=1
            addTris = [addTris;[patientEdge(1:roww(1),:),repmat(rootEdges(1,1),roww(1),1)]];
        end

    elseif  roww(end) < length(patientEdge)/2 && roww(1)~=1
       s = find(roww > length(patientEdge)/2);
       endIdx = s(end);        % roww�����һ������С��length(patientEdge)/2��������
       addTris = [addTris;[patientEdge(roww(endIdx):length(patientEdge),:),repmat(rootEdges(endIdx,2),length(patientEdge)-roww(endIdx)+1,1)]];
       addTris = [addTris;[patientEdge(1:roww(end),:),repmat(rootEdges(endIdx,2),roww(end),1)]]; 
       addTris = [addTris;[patientEdge(roww(end):roww(1)-1,:),repmat(rootEdges(1,1),roww(1)-roww(end),1)]]; 
    end

else       %12hit
    s = find(roww < length(patientEdge)/2);      
    start = s(1);            % roww�е�һ������С��length(patientEdge)/2��������
   
   for j = start:length(roww)-1
        if roww(j) ~= roww(j+1) 
            temp11 = roww(j):(roww(j+1)-1);
            temp1 = patientEdge(temp11,:);
            temp2 = repmat(rootEdges(j,2),roww(j+1)-roww(j),1);
            temp = [temp1, temp2];
            addTris = [addTris; temp];   
        end
   end 

   if start>2     % 12hit
        for j = 1:start-2           %1~4
             if roww(j) ~= roww(j+1) 
                 temp11 = roww(j):roww(j+1)-1;
                 temp1 = patientEdge(temp11,:);
                 temp2 = repmat(rootEdges(j,2),roww(j+1)-roww(j),1);
                 temp = [temp1,temp2];
                addTris = [addTris;temp];   
             end
        end
   end

   if roww(end)<roww(1) % 12hit
      addTris = [addTris;[patientEdge(roww(end):roww(1)-1,:),repmat(rootEdges(end,2),roww(1)-roww(end),1)]];
   end

   if roww(start-1)<=length(patientEdge)    %12hit
         addTris = [addTris;[patientEdge(roww(start-1):length(patientEdge),:),repmat(rootEdges(start-1,2),length(patientEdge)-roww(start-1)+1,1)]];
   end

   addTris = [addTris;[patientEdge(1:roww(start)-1,:),repmat(rootEdges(start,1),roww(start)-1,1)]];

end

newTris = [allTris; addTris];

DatWriteIdxs('addTris.dat', addTris);
writeOBJ('����ǰ�ĺϲ�����.obj', allVers, newTris);
OBJwriteVertices('centerPatient.obj', centerPatient);

%% 5  ���ο��Ʋ���
index_crown_change = find(patientTooth.vertex(:,2)>(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-1 ...
                        & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+2);
patientTransform = patientTooth.vertex(index_crown_change,:);

%���������ֱ��γ�����
mergedToothVers = allVers;
mergeRegionIdx =  find(mergedToothVers(:,2)>(centerPatient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-1 ...
    &mergedToothVers(:,2)<=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+1);
mergeRegionVers = mergedToothVers(mergeRegionIdx,:);

% �������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�movedRootTooth2����crownforchange��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
indices = 1:size(mergedToothVers,1);

% �ϲ������У�����Ҫ������εĶ������������һ����������
exterior = indices(mergedToothVers(:,2)<(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-0.6...
    |mergedToothVers(:,2)>=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+3);
 
OBJwriteVertices('mergeRegionVers.obj', mergeRegionVers);
OBJwriteVertices('patientTransform.obj', patientTransform);
exteriorVers = mergedToothVers(exterior, :);
OBJwriteVertices('exteriorVers.obj', exteriorVers);

DatWriteIdxs('mergeRegionIdx.dat', mergeRegionIdx);
DatWriteIdxs('exterior.dat', exterior');

    end
   
else
    fdi = textread('FDILower__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
        disp( '�ò���û�д����� ');
        return;
    else

%%  1.�и������
        namestr1 = ['toothLower_',num2str(toothIdx-1),'.','obj'];
        patientTooth = Read_Obj(namestr1);
        namestr2 = ['gumlineLower_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);
        axis = ReadObj('AXISLower_.obj');
        axisPatient = axis(3*(toothIdx-1)+1:3*toothIdx,:);
        
        axisStandardVers = [axisStandard.x; axisStandard.y; axisStandard.z];
OBJwriteVertices('axisPatient.obj', axisPatient);
OBJwriteVertices('axisStandard.obj', axisStandardVers);
OBJwriteVertices('gumline.obj', gumline);
writeOBJ('patientTooth.obj',patientTooth.vertex, patientTooth.face);
writeOBJ('crownTooth.obj',crownTooth.vertex, crownTooth.face);
writeOBJ('rootTooth.obj',rootTooth.vertices, rootTooth.faces);

         %����
        centerPatient = mean(patientTooth.vertex);y = min(gumline(:,2));
        p_patient =[centerPatient(1),y(1),centerPatient(3)];patientYdir = axisPatient(2,:);patientTooth.vertex = patientTooth.vertex;patientTooth.face = patientTooth.face;
        cutPatientIdx = find(patientTooth.vertex(:,2)>=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - ... 
            p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+0.5);
        cutPatientCrownVers = patientTooth.vertex(cutPatientIdx,:);
         a1 = ismember(patientTooth.face(:,1),cutPatientIdx);
        aa1 = find(a1==1);
        a2 = ismember(patientTooth.face(:,2),cutPatientIdx);
        aa2 = find(a2==1);
        a3 = ismember(patientTooth.face(:,3),cutPatientIdx);
        aa3 = find(a3==1);
        [c1,~,~] = intersect(aa1,aa2);
        triIdxInCutPatient = intersect(c1,aa3);
        trisInCutPatient = patientTooth.face(triIdxInCutPatient,:);
        
        %�ұ߽�
        [raw_edges_list] = query_edges_list(patientTooth.face(triIdxInCutPatient,:),'sorted');
        [~,~,iu] = unique(sort(raw_edges_list,2),'rows');
        [i3,i4] = histc(iu,unique(iu));
        lone_edges_idx_vect = i3(i4) == 1;
        bdryEdges = unique(raw_edges_list(lone_edges_idx_vect,:),'rows');
 
        %�߽�㣺
        edgeVerIdx = unique([bdryEdges(:,1);bdryEdges(:,2)]);
        
        %�߽�������ڵ��е�λ��
        [~,tempRowIdx,~] = intersect(cutPatientIdx,edgeVerIdx);
        edge_p_b = cutPatientCrownVers(tempRowIdx,:);
        patientBdryEdges_newRep = zeros(size(bdryEdges));
        for i = 1:length(tempRowIdx)
            index = find(bdryEdges == edgeVerIdx(i));
            patientBdryEdges_newRep(index) = tempRowIdx(i);
        end
        patientEdge = sotr_edge(patientBdryEdges_newRep,1);
        
      
        patientCutTris = ones (size(trisInCutPatient));     % ��������ʾ�Ĳ����и����������Ƭ��
        for j = 1:length(cutPatientCrownVers)
            index = find( trisInCutPatient == cutPatientIdx(j));
            patientCutTris(sub2ind(size(patientCutTris), index)) = j;       
        end

        writeOBJ('�и��Ĳ�����������.obj', cutPatientCrownVers, patientCutTris)

        
%%  2. ���롪�����ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж���

%  2.1 ȷ����һ����תƽ���������     
        centerCrown = mean(crownTooth.vertex);
        rootYdir = axisStandard.y;
        tempCrownIdx = find(crownTooth.vertex(:,2)>=(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(:,1) ...
                    - centerCrown(1))+rootYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/rootYdir(2)));
       
        for i =1:length(tempCrownIdx)
            temp11 = crownTooth.vertex(tempCrownIdx(i), :);
            temp1 = temp11 - centerCrown;
            crownValue(i) = abs( dot(temp1, rootYdir));
        end
      
        maxIdx = find(crownValue == max(crownValue));
        point_crown =  crownTooth.vertex(tempCrownIdx(maxIdx),:);


        %��׼����
        centerRoot = mean(rootTooth.vertices);
        tempRootIdx = find(rootTooth.vertices(:,2)>=(centerRoot(2) - (rootYdir(1)*(rootTooth.vertices(:,1) ...
                         - centerRoot(1))+rootYdir(3)*(rootTooth.vertices(:,3) - centerRoot(3)))/rootYdir(2)));
        
        for i =1:length(tempRootIdx)
            temp11 = rootTooth.vertices(tempRootIdx(i),:);
            temp1 = temp11-centerRoot;
            rootValue(i) = abs(dot(temp1,rootYdir));
        end
        
        maxIdx = find(rootValue == max(rootValue));
        valueRootMaxIdx = tempRootIdx(maxIdx);
        point_root = rootTooth.vertices(valueRootMaxIdx,:);

        
%       2.2 ��һ����תƽ��
        movingVec = point_crown - point_root;
        movingMat = repmat(movingVec,length(rootTooth.vertices),1);
        movedRootTooth0 = rootTooth.vertices + movingMat;
        
        temp = inv([axisStandard.x;axisStandard.y;axisStandard.z]);
        R = temp * axisPatient;
        movingVec = centerPatient - centerCrown*R;
        movingMat = repmat(movingVec,length(rootTooth.vertices),1);
        movedRootTooth1 = movedRootTooth0 *R + movingMat;
 
       OBJwriteVertices('��һ����תƽ�ƺ�Ĵ�����׼��.obj', movedRootTooth1);
 
        %����仯���ݵ���̬
        index_1 = find(patientTooth.vertex(:,2)<=(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2)) ...
                & patientTooth.vertex(:,2)>=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2)+0.5));
        point1 = patientTooth.vertex(index_1,:);
        index_2 = find(movedRootTooth1(:,2)<=(centerPatient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2))  ...
                & movedRootTooth1(:,2)>=(p_patient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2)+0.5));
        point2 = movedRootTooth1(index_2,:);
 
%       2.4 icp����
        for i = 1:length(point2)
           [minValue,startEdgeIdx]=mindis(point1,point2(i,:),1);
           row(i) = startEdgeIdx;
        end
        point11 = point1(row,:);
        [R,t,~,~,~,~] = icp(point11,point2);
        
%       2.5 �ڶ�����תƽ��
        movedRootTooth2 = bsxfun(@plus,movedRootTooth1,t);
        
        save('movedRootTooth2.mat', 'movedRootTooth2');
        OBJwriteVertices('�ڶ�����תƽ�ƺ�Ĵ�����׼��.obj', movedRootTooth2);


%% 3. �и������׼��  

%       3.1 ȷ���и�㡣
rootCutIdx =  find(movedRootTooth2(:,2)<(p_patient(2) - ( patientYdir(1)*(movedRootTooth2(:,1) ...
                 - p_patient(1))+ patientYdir(3)*(movedRootTooth2(:,3) - p_patient(3)))/ patientYdir(2))+0.5);%�г���������������ֵ%ע�⣡���ղ������ڷ�����
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

writeOBJ('������Ƭǰ�ĺϲ�����.obj', allVers, allTris);
OBJwriteVertices('�и���������.obj', rootCutVers);
writeOBJ('�и���������.obj', rootCutVers, trisInCutRoot);


%% 4. �ϲ���������Ƭ

%5.20����
%1.���еı߽�㰴��һ����������
re1 = patientEdge(1,1);                             % �������ڱ�Ե���һ�����������
rp11 = allVers(re1,:);
rbeVersIdx = rootBdryEdges_newRep(:,1);                         % ������Ե�������
[~,startEdgeIdx] = mindis(allVers(rbeVersIdx,:), rp11, 1);       % ������Ե���о������ڱ�Ե��һ����������������

% �ж�[patientEdge(1,��),Edge_t(r,1)]de ������[Edge_t(r,��),patientEdge(1,1)]�ķ����Ƿ���ͬ
% Ϊ�˷�ֹ���������������Edge_t(r,1)����ĵڶ�����
afterStartIdx = rootBdryEdges_newRep(startEdgeIdx,2);       % ��ʼ�ߵĺ󶥵��������
[row1,col1] = ind2sub(size(rootBdryEdges_newRep),find(rootBdryEdges_newRep == afterStartIdx));     

rr = row1(ismember(row1,startEdgeIdx)==0);      % row1����Ԫ���в�����startIdx�е��Ǹ�
rm = col1(ismember(row1,startEdgeIdx)==0);      % 1����2������ѡȡ��������Ӧ��λ�á�

pe3 = patientEdge(3,1);
pe1 = patientEdge(1,1);
pp1 = allVers(pe1,:);
pp3 = allVers(pe3,:) ;
rp1 = allVers(rootBdryEdges_newRep(startEdgeIdx,1),:);
rp3 = allVers(rootBdryEdges_newRep(rr,3-rm),:);
rCenter = mean(rootEdgeVers);

tri11 = pp3 - pp1;
tri12 = rp1 - pp3;
tri21 = rp3 - rp1;
tri22 = pp1 - rp3;

s_b = dot(pp1 - rCenter, cross(tri11,tri12));
s_t = dot(rp1 - rCenter, cross(tri21,tri22));

if s_b * s_t <= 0   % 11hit  ��ͬ�� 
    rootEdges = sotr_edge(rootBdryEdges_newRep, startEdgeIdx);
else            % 12hit
    E_t = [rootBdryEdges_newRep(:,2),rootBdryEdges_newRep(:,1)];
    rootEdges = sotr_edge(E_t,startEdgeIdx);
end


%   4.2 
addTris = [];
middle = [];


for j = 1:length(rootEdges)
    pe1 = allVers(rootEdges(j,1),:);
    p2 = allVers(rootEdges(j,2),:);
    middle(j, :) = (pe1+p2)/2;
    [~,ro]=mindis(allVers(patientEdge(:,1),:),(pe1+p2)/2,1);
    roww(j) = ro;
    addTris = [addTris;[rootEdges(j,:),patientEdge(ro,1)]];
end

OBJwriteVertices('middle.obj',middle);


%   4.3. 
if roww(1) < length(patientEdge)/2    % 11hit

    for j = 1:length(roww)-1
        if roww(j) ~= roww(j+1) 
            addTris = [addTris;[patientEdge(roww(j):roww(j+1)-1,:),repmat(rootEdges(j,2),roww(j+1)-roww(j),1)]];
        end
    end

    if roww(end) == length(patientEdge)
        addTris = [addTris; [patientEdge(roww(end),:), rootEdges(1,1)]];

        if roww(1)~=1
            addTris = [addTris;[patientEdge(1:roww(1),:),repmat(rootEdges(1,1),roww(1),1)]];
        end

    elseif roww(end) > length(patientEdge)/2 && roww(end)<length(patientEdge)
        addTris = [addTris;[patientEdge(roww(end):length(patientEdge),:),repmat(rootEdges(1,1),length(patientEdge)-roww(end)+1,1)]];        
        if roww(1)~=1
            addTris = [addTris;[patientEdge(1:roww(1),:),repmat(rootEdges(1,1),roww(1),1)]];
        end

    elseif  roww(end) < length(patientEdge)/2 && roww(1)~=1
       s = find(roww > length(patientEdge)/2);
       endIdx = s(end);        % roww�����һ������С��length(patientEdge)/2��������
       addTris = [addTris;[patientEdge(roww(endIdx):length(patientEdge),:),repmat(rootEdges(endIdx,2),length(patientEdge)-roww(endIdx)+1,1)]];
       addTris = [addTris;[patientEdge(1:roww(end),:),repmat(rootEdges(endIdx,2),roww(end),1)]]; 
       addTris = [addTris;[patientEdge(roww(end):roww(1)-1,:),repmat(rootEdges(1,1),roww(1)-roww(end),1)]]; 
    end

else       %12hit
    s = find(roww < length(patientEdge)/2);      
    start = s(1);            % roww�е�һ������С��length(patientEdge)/2��������
   
   for j = start:length(roww)-1
        if roww(j) ~= roww(j+1) 
            temp11 = roww(j):(roww(j+1)-1);
            temp1 = patientEdge(temp11,:);
            temp2 = repmat(rootEdges(j,2),roww(j+1)-roww(j),1);
            temp = [temp1, temp2];
            addTris = [addTris; temp];   
        end
   end 

   if start>2     % 12hit
        for j = 1:start-2           %1~4
             if roww(j) ~= roww(j+1) 
                 temp11 = roww(j):roww(j+1)-1;
                 temp1 = patientEdge(temp11,:);
                 temp2 = repmat(rootEdges(j,2),roww(j+1)-roww(j),1);
                 temp = [temp1,temp2];
                addTris = [addTris;temp];   
             end
        end
   end

   if roww(end)<roww(1) % 12hit
      addTris = [addTris;[patientEdge(roww(end):roww(1)-1,:),repmat(rootEdges(end,2),roww(1)-roww(end),1)]];
   end

   if roww(start-1)<=length(patientEdge)    %12hit
         addTris = [addTris;[patientEdge(roww(start-1):length(patientEdge),:),repmat(rootEdges(start-1,2),length(patientEdge)-roww(start-1)+1,1)]];
   end

   addTris = [addTris;[patientEdge(1:roww(start)-1,:),repmat(rootEdges(start,1),roww(start)-1,1)]];

end

newTris = [allTris; addTris];

writeOBJ('����ǰ�ĺϲ�����.obj', allVers, newTris);
OBJwriteVertices('centerPatient.obj', centerPatient);


%% 5  ���ο��Ʋ���
index_crown_change = find(patientTooth.vertex(:,2)<(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+1 ...
                        & patientTooth.vertex(:,2)>=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-2);
patientTransform = patientTooth.vertex(index_crown_change,:);

%���������ֱ��γ�����
mergedToothVers = allVers;
mergeRegionIdx =  find(mergedToothVers(:,2)<(centerPatient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+1 ...
    &mergedToothVers(:,2)>=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-1);
mergeRegionVers = mergedToothVers(mergeRegionIdx,:);

% �������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�movedRootTooth2����crownforchange��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
indices = 1:size(mergedToothVers,1);

% �ϲ������У�����Ҫ������εĶ������������һ����������
exterior = indices(mergedToothVers(:,2)>(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+1.5...
    |mergedToothVers(:,2)<=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-3);
 
OBJwriteVertices('mergeRegionVers.obj', mergeRegionVers);
OBJwriteVertices('patientTransform.obj', patientTransform);
exteriorVers = mergedToothVers(exterior, :);
OBJwriteVertices('exteriorVers.obj', exteriorVers);

DatWriteIdxs('mergeRegionIdx.dat', mergeRegionIdx);
DatWriteIdxs('exterior.dat', exterior');
 

    end
    

        
     
end
 

%% 6. ����
[omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(mergedToothVers,1), newTris, exterior);

%[bi_L,bi_U,bi_P,bi_Q,bi_D,bi_S,bi_M] = biharm_factor_system(mergedToothVers,newTris, 'ext','voronoi', 'no_flatten',omega,N0,N1);
[A, bi_S] =  biharm_factor_system_modified(mergedToothVers,newTris,omega,N0,N1);

%�ҵ��������ָ����ڱ��β�������ĵ㣬�������ϵĵ��������ڴ�
finalVers = mergedToothVers;
for i = 1:length(mergeRegionVers)
   [minValue,startEdgeIdx] = mindis(patientTransform, mergeRegionVers(i,:),1);
   minvalue(i) = minValue;  row(i) = startEdgeIdx;
   finalVers(mergeRegionIdx(i),:) = patientTransform(startEdgeIdx,:);
end

 
BZ1 = zeros(size(mergedToothVers,1),3);
finalVers = biharm_solve_with_factor_modified(A, bi_S, finalVers, omega, N0, N1);
 
writeOBJ('��������.obj',finalVers, newTris)


%% 7. ��������Ƭ����
focusVers = finalVers(mergeRegionIdx, :);
OBJwriteVertices('focusVers.obj', focusVers);

focusCenter = mean(focusVers);      % ����ƴ�Ӵ��������ж�����е㡣
deltaHeight = 3;


disp('finished.');
