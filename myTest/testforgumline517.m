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
addpath([functiondir '牙轴对齐（标准牙和病人牙冠对齐）'])
addpath([functiondir 'siofmodeltooth'])
addpath([functiondir 'attachments'])
addpath([functiondir 'modelStadard'])
addpath([functiondir 'testdata'])
addpath([functiondir 'mesh process'])
addpath([functiondir 'testdata(4)'])
% % 加载标准牙
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
crownTooth = crown.model;          % 标准牙冠
        
        
s = fix(x/10);  %取整
g = mod(x,10);%取余

    
if s == 1 || s == 2
    fdi = textread('FDIUpper__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
        disp( '该病例没有此牙齿 ');
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

%%  1.切割病人牙冠

%      1.2 确定病人牙冠切割顶点。
        centerPatient = mean(patientTooth.vertex);
        y = max(gumline(:,2));
        p_patient =[centerPatient(1),y(1),centerPatient(3)];
        patientYdir = axisPatient(2,:);
 
        cutPatientIdx = find(patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-2);
        cutPatientCrownVers = patientTooth.vertex(cutPatientIdx,:);
        
%       1.3 确定病人切割牙冠边界、三角片。
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
        
        %边界点：
        edgeVerIdx = unique([bdryEdges(:,1);bdryEdges(:,2)]);
        
        %边界点在牙冠点中的位置
        [~,tempRowIdx,~] = intersect(cutPatientIdx,edgeVerIdx);
        
 %           老索引——原病人网格中的点索引； 新索引——切割之后的病人网格中的点索引；
        edge_p_b = cutPatientCrownVers(tempRowIdx,:);
        patientBdryEdges_newRep = zeros(size(bdryEdges));
        for i = 1:length(tempRowIdx)
            index = find(bdryEdges == edgeVerIdx(i));
            patientBdryEdges_newRep(index) = tempRowIdx(i);
        end
        
        patientEdge = sotr_edge(patientBdryEdges_newRep,1);    
        patientCutTris = ones (size(trisInCutPatient)); % 新索引表示的病人切割网格的三角片。
        
        for j = 1:length(cutPatientCrownVers)
            index = find( trisInCutPatient == cutPatientIdx(j));
            patientCutTris(sub2ind(size(patientCutTris), index)) = j;
        end
        
        
writeOBJ('切割后的病人牙齿网格.obj', cutPatientCrownVers, patientCutTris)



%% 2. 对齐——利用标准牙和病人牙齿的中心点以及三轴进行对齐

%  2.1 确定第一次旋转平移所需参数
        centerCrown = mean(crownTooth.vertex);
        rootYdir = axisStandard.y;
        tempCrownIdx = find(crownTooth.vertex(:,2)<=(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(:,1) ...
                    - centerCrown(1))+rootYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/rootYdir(2)));
        crownValue = abs(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(tempCrownIdx,1)-centerCrown(1))...
                     +rootYdir(3)*(crownTooth.vertex(tempCrownIdx,3)-centerCrown(3)))/rootYdir(2));
        maxIdx = find(crownValue == max(crownValue));
        point_crown =  crownTooth.vertex(tempCrownIdx(maxIdx),:);

        %标准牙根
        centerRoot = mean(rootTooth.vertices);
        tempRootIdx = find(rootTooth.vertices(:,2)<=(centerRoot(2) - (rootYdir(1)*(rootTooth.vertices(:,1) ...
                         - centerRoot(1))+rootYdir(3)*(rootTooth.vertices(:,3) - centerRoot(3)))/rootYdir(2)));
        rootValue = abs(centerRoot(2) - (rootYdir(1)*(rootTooth.vertices(tempRootIdx,1)-centerRoot(1))...
                   +rootYdir(3)*(rootTooth.vertices(tempRootIdx,3)-centerRoot(3)))/rootYdir(2));
        maxIdx = find(rootValue == max(rootValue));
        point_tooth = rootTooth.vertices(tempRootIdx(maxIdx),:);
        centcrownintooth = centerCrown+point_tooth-point_crown;%切后标准根的边缘中心
        
        
        %%利用标准牙和病人牙齿的中心点以及三轴进行对齐，对齐后按照病人牙冠的切割位置进行切割
        temp = inv([axisStandard.x;axisStandard.y;axisStandard.z]);
        R = temp * axisPatient;
        movingVec = centerPatient - centcrownintooth*R;
        movingMat = repmat(movingVec,length(rootTooth.vertices),1);
        movedRootTooth1 = (rootTooth.vertices) *R + movingMat;
        
        OBJwriteVertices('第一次旋转平移后的带根标准牙.obj', movedRootTooth1);


        %整体变化牙齿的形态
        index_1 = find(patientTooth.vertex(:,2)>=(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))...
                & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+0.5);
        point1 = patientTooth.vertex(index_1,:);
        
        index_2 = find(movedRootTooth1(:,2)>=(centerPatient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2))...
                & movedRootTooth1(:,2)<=(p_patient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2))+0.5);
        point2 = movedRootTooth1(index_2,:);
        
        
%       2.4 icp计算
        for i = 1:length(point2)
           [minValue,startEdgeIdx]=mindis(point1,point2(i,:),1); 
           row(i) = startEdgeIdx;
        end
        
        point11 = point1(row,:);
        [R,t,~,~,~,~] = icp(point11,point2);
        
%       2.5 第二次旋转平移
        movedRootTooth2 = bsxfun(@plus,movedRootTooth1*R,t);

        
        save('movedRootTooth2.mat', 'movedRootTooth2');
        OBJwriteVertices('第二次旋转平移后的带根标准牙.obj', movedRootTooth2);

        

%% 3. 切割带根标准牙  

%       3.1 确定切割顶点。
rootCutIdx =  find(movedRootTooth2(:,2)>(p_patient(2) - ( patientYdir(1)*(movedRootTooth2(:,1) ...
                 - p_patient(1))+ patientYdir(3)*(movedRootTooth2(:,3) - p_patient(3)))/ patientYdir(2))-2);%切出来的牙根的索引值%注意！按照病人牙冠方向切
rootCutVers = movedRootTooth2(rootCutIdx ,:);

a_t1 = ismember(rootTooth.faces(:,1),rootCutIdx);
aa_t1 = find(a_t1==1);
a_t2 = ismember(rootTooth.faces(:,2),rootCutIdx);
aa_t2 = find(a_t2==1);
a_t3 = ismember(rootTooth.faces(:,3),rootCutIdx);
aa_t3 = find(a_t3==1);
[c_t1,~,~] = intersect(aa_t1,aa_t2);
triIdxInCutRoot = intersect(c_t1,aa_t3);     % 列向量，root中三个顶点索引都能在cutRoot中找到的三角片的索引。

%       3.2 确定切割网格的边界和三角片。
[raw_edges_list_t] = query_edges_list(rootTooth.faces(triIdxInCutRoot,:),'sorted');
[~,~,iu_t] = unique(sort(raw_edges_list_t,2),'rows');
[i3,i4] = histc(iu_t,unique(iu_t));
lone_edges_idx_vect_t = i3(i4) == 1;
rootBdryEdges = unique(raw_edges_list_t(lone_edges_idx_vect_t,:),'rows');       %  用老索引表示的切割牙根边界边。

trisInCutRoot = rootTooth.faces(triIdxInCutRoot,:);


%       3.3 确定切割网格三角片、边界映射到合并网格中的三角片、边界。
rootCutTrisInMergeMesh = ones (size(trisInCutRoot));           % 新索引表示的cutRoot中的三角片。
for j = 1:length(rootCutVers)
    index = find( trisInCutRoot == rootCutIdx(j));
    rootCutTrisInMergeMesh(sub2ind(size(rootCutTrisInMergeMesh), index)) = j+length(cutPatientCrownVers);
end

%边界点：
rootEdgeVersIdx = unique([rootBdryEdges(:,1);rootBdryEdges(:,2)]);  % 点索引向量，边界中所有点的索引，无重复

%边界点在牙根点中的位置
[~,tempRowIdx,~] = intersect(rootCutIdx,rootEdgeVersIdx);

rootEdgeVers = rootCutVers(tempRowIdx,:);
tempRowIdx = tempRowIdx + length(cutPatientCrownVers);

rootBdryEdges_newRep = zeros(size(rootBdryEdges));          % 用新索引表示的边界边。
for i = 1:length(tempRowIdx)
    index = find(rootBdryEdges == rootEdgeVersIdx(i));
    rootBdryEdges_newRep(index) = tempRowIdx(i);
end


allTris = [patientCutTris; rootCutTrisInMergeMesh];            % 合并网格的顶点。
allVers = [cutPatientCrownVers; rootCutVers];       % 合并网格三角片，未整理的


save('rootBdryEdges_newRep.mat', 'rootBdryEdges_newRep');
save('patientEdge.mat', 'patientEdge');     
save('allVers.mat', 'allVers');
save('allTris.mat', 'allTris');
save('rootEdgeVers.mat', 'rootEdgeVers');

writeOBJ('补三角片前的合并网格.obj', allVers, allTris);
OBJwriteVertices('切割牙根顶点.obj', rootCutVers);
writeOBJ('切割牙根网格.obj', rootCutVers, trisInCutRoot);


%% 4. 合并网格补三角片

%5.20更新
%1.所有的边界点按照一个方向排序
re1 = patientEdge(1,1);                             % 病人牙冠边缘点第一个点的索引。
rp11 = allVers(re1,:);
rbeVersIdx = rootBdryEdges_newRep(:,1);                         % 牙根边缘点的索引
[~,startEdgeIdx] = mindis(allVers(rbeVersIdx,:), rp11, 1);       % 牙根边缘点中距离牙冠边缘第一个点最近点的索引。

% 判断[patientEdge(1,：),Edge_t(r,1)]de 方向与[Edge_t(r,：),patientEdge(1,1)]的方向是否相同
% 为了防止个别特殊情况，找Edge_t(r,1)后面的第二个点
afterStartIdx = rootBdryEdges_newRep(startEdgeIdx,2);       % 起始边的后顶点的索引。
[row1,col1] = ind2sub(size(rootBdryEdges_newRep),find(rootBdryEdges_newRep == afterStartIdx));     

rr = row1(ismember(row1,startEdgeIdx)==0);      % row1两个元素中不等于startIdx中的那个
rm = col1(ismember(row1,startEdgeIdx)==0);      % 1或者2，上面选取的索引对应的位置。

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

if s_b * s_t <= 0   % 11hit  不同号 
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
       endIdx = s(end);        % roww中最后一个满足小于length(patientEdge)/2的索引；
       addTris = [addTris;[patientEdge(roww(endIdx):length(patientEdge),:),repmat(rootEdges(endIdx,2),length(patientEdge)-roww(endIdx)+1,1)]];
       addTris = [addTris;[patientEdge(1:roww(end),:),repmat(rootEdges(endIdx,2),roww(end),1)]]; 
       addTris = [addTris;[patientEdge(roww(end):roww(1)-1,:),repmat(rootEdges(1,1),roww(1)-roww(end),1)]]; 
    end

else       %12hit
    s = find(roww < length(patientEdge)/2);      
    start = s(1);            % roww中第一个满足小于length(patientEdge)/2的索引；
   
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
writeOBJ('变形前的合并网格.obj', allVers, newTris);
OBJwriteVertices('centerPatient.obj', centerPatient);

%% 5  变形控制部分
index_crown_change = find(patientTooth.vertex(:,2)>(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-1 ...
                        & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+2);
patientTransform = patientTooth.vertex(index_crown_change,:);

%牙冠牙根分别形成网格
mergedToothVers = allVers;
mergeRegionIdx =  find(mergedToothVers(:,2)>(centerPatient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-1 ...
    &mergedToothVers(:,2)<=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+1);
mergeRegionVers = mergedToothVers(mergeRegionIdx,:);

% 病人牙冠部分多截一部分出来，用以确定标准牙根变形后点的位置，找到movedRootTooth2中离crownforchange中最近点的位置，并将这些点变为新的点，即牙冠点的位置
indices = 1:size(mergedToothVers,1);

% 合并网格中，不需要参与变形的顶点的索引，是一个行向量。
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
        disp( '该病例没有此牙齿 ');
        return;
    else

%%  1.切割病人牙冠
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

         %处理
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
        
        %找边界
        [raw_edges_list] = query_edges_list(patientTooth.face(triIdxInCutPatient,:),'sorted');
        [~,~,iu] = unique(sort(raw_edges_list,2),'rows');
        [i3,i4] = histc(iu,unique(iu));
        lone_edges_idx_vect = i3(i4) == 1;
        bdryEdges = unique(raw_edges_list(lone_edges_idx_vect,:),'rows');
 
        %边界点：
        edgeVerIdx = unique([bdryEdges(:,1);bdryEdges(:,2)]);
        
        %边界点在牙冠点中的位置
        [~,tempRowIdx,~] = intersect(cutPatientIdx,edgeVerIdx);
        edge_p_b = cutPatientCrownVers(tempRowIdx,:);
        patientBdryEdges_newRep = zeros(size(bdryEdges));
        for i = 1:length(tempRowIdx)
            index = find(bdryEdges == edgeVerIdx(i));
            patientBdryEdges_newRep(index) = tempRowIdx(i);
        end
        patientEdge = sotr_edge(patientBdryEdges_newRep,1);
        
      
        patientCutTris = ones (size(trisInCutPatient));     % 新索引表示的病人切割网格的三角片。
        for j = 1:length(cutPatientCrownVers)
            index = find( trisInCutPatient == cutPatientIdx(j));
            patientCutTris(sub2ind(size(patientCutTris), index)) = j;       
        end

        writeOBJ('切割后的病人牙齿网格.obj', cutPatientCrownVers, patientCutTris)

        
%%  2. 对齐——利用标准牙和病人牙齿的中心点以及三轴进行对齐

%  2.1 确定第一次旋转平移所需参数     
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


        %标准牙根
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

        
%       2.2 第一次旋转平移
        movingVec = point_crown - point_root;
        movingMat = repmat(movingVec,length(rootTooth.vertices),1);
        movedRootTooth0 = rootTooth.vertices + movingMat;
        
        temp = inv([axisStandard.x;axisStandard.y;axisStandard.z]);
        R = temp * axisPatient;
        movingVec = centerPatient - centerCrown*R;
        movingMat = repmat(movingVec,length(rootTooth.vertices),1);
        movedRootTooth1 = movedRootTooth0 *R + movingMat;
 
       OBJwriteVertices('第一次旋转平移后的带根标准牙.obj', movedRootTooth1);
 
        %整体变化牙齿的形态
        index_1 = find(patientTooth.vertex(:,2)<=(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2)) ...
                & patientTooth.vertex(:,2)>=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2)+0.5));
        point1 = patientTooth.vertex(index_1,:);
        index_2 = find(movedRootTooth1(:,2)<=(centerPatient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2))  ...
                & movedRootTooth1(:,2)>=(p_patient(2) - (patientYdir(1)*(movedRootTooth1(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth1(:,3) - p_patient(3)))/patientYdir(2)+0.5));
        point2 = movedRootTooth1(index_2,:);
 
%       2.4 icp计算
        for i = 1:length(point2)
           [minValue,startEdgeIdx]=mindis(point1,point2(i,:),1);
           row(i) = startEdgeIdx;
        end
        point11 = point1(row,:);
        [R,t,~,~,~,~] = icp(point11,point2);
        
%       2.5 第二次旋转平移
        movedRootTooth2 = bsxfun(@plus,movedRootTooth1,t);
        
        save('movedRootTooth2.mat', 'movedRootTooth2');
        OBJwriteVertices('第二次旋转平移后的带根标准牙.obj', movedRootTooth2);


%% 3. 切割带根标准牙  

%       3.1 确定切割顶点。
rootCutIdx =  find(movedRootTooth2(:,2)<(p_patient(2) - ( patientYdir(1)*(movedRootTooth2(:,1) ...
                 - p_patient(1))+ patientYdir(3)*(movedRootTooth2(:,3) - p_patient(3)))/ patientYdir(2))+0.5);%切出来的牙根的索引值%注意！按照病人牙冠方向切
rootCutVers = movedRootTooth2(rootCutIdx ,:);

a_t1 = ismember(rootTooth.faces(:,1),rootCutIdx);
aa_t1 = find(a_t1==1);
a_t2 = ismember(rootTooth.faces(:,2),rootCutIdx);
aa_t2 = find(a_t2==1);
a_t3 = ismember(rootTooth.faces(:,3),rootCutIdx);
aa_t3 = find(a_t3==1);
[c_t1,~,~] = intersect(aa_t1,aa_t2);
triIdxInCutRoot = intersect(c_t1,aa_t3);     % 列向量，root中三个顶点索引都能在cutRoot中找到的三角片的索引。

%       3.2 确定切割网格的边界和三角片。
[raw_edges_list_t] = query_edges_list(rootTooth.faces(triIdxInCutRoot,:),'sorted');
[~,~,iu_t] = unique(sort(raw_edges_list_t,2),'rows');
[i3,i4] = histc(iu_t,unique(iu_t));
lone_edges_idx_vect_t = i3(i4) == 1;
rootBdryEdges = unique(raw_edges_list_t(lone_edges_idx_vect_t,:),'rows');       %  用老索引表示的切割牙根边界边。

trisInCutRoot = rootTooth.faces(triIdxInCutRoot,:);


%       3.3 确定切割网格三角片、边界映射到合并网格中的三角片、边界。
rootCutTrisInMergeMesh = ones (size(trisInCutRoot));           % 新索引表示的cutRoot中的三角片。
for j = 1:length(rootCutVers)
    index = find( trisInCutRoot == rootCutIdx(j));
    rootCutTrisInMergeMesh(sub2ind(size(rootCutTrisInMergeMesh), index)) = j+length(cutPatientCrownVers);
end


%边界点：
rootEdgeVersIdx = unique([rootBdryEdges(:,1);rootBdryEdges(:,2)]);  % 点索引向量，边界中所有点的索引，无重复

%边界点在牙根点中的位置
[~,tempRowIdx,~] = intersect(rootCutIdx,rootEdgeVersIdx);

rootEdgeVers = rootCutVers(tempRowIdx,:);
tempRowIdx = tempRowIdx + length(cutPatientCrownVers);

rootBdryEdges_newRep = zeros(size(rootBdryEdges));          % 用新索引表示的边界边。
for i = 1:length(tempRowIdx)
    index = find(rootBdryEdges == rootEdgeVersIdx(i));
    rootBdryEdges_newRep(index) = tempRowIdx(i);
end


allTris = [patientCutTris; rootCutTrisInMergeMesh];            % 合并网格的顶点。
allVers = [cutPatientCrownVers; rootCutVers];       % 合并网格三角片，未整理的


save('rootBdryEdges_newRep.mat', 'rootBdryEdges_newRep');
save('patientEdge.mat', 'patientEdge');     
save('allVers.mat', 'allVers');
save('allTris.mat', 'allTris');
save('rootEdgeVers.mat', 'rootEdgeVers');

writeOBJ('补三角片前的合并网格.obj', allVers, allTris);
OBJwriteVertices('切割牙根顶点.obj', rootCutVers);
writeOBJ('切割牙根网格.obj', rootCutVers, trisInCutRoot);


%% 4. 合并网格补三角片

%5.20更新
%1.所有的边界点按照一个方向排序
re1 = patientEdge(1,1);                             % 病人牙冠边缘点第一个点的索引。
rp11 = allVers(re1,:);
rbeVersIdx = rootBdryEdges_newRep(:,1);                         % 牙根边缘点的索引
[~,startEdgeIdx] = mindis(allVers(rbeVersIdx,:), rp11, 1);       % 牙根边缘点中距离牙冠边缘第一个点最近点的索引。

% 判断[patientEdge(1,：),Edge_t(r,1)]de 方向与[Edge_t(r,：),patientEdge(1,1)]的方向是否相同
% 为了防止个别特殊情况，找Edge_t(r,1)后面的第二个点
afterStartIdx = rootBdryEdges_newRep(startEdgeIdx,2);       % 起始边的后顶点的索引。
[row1,col1] = ind2sub(size(rootBdryEdges_newRep),find(rootBdryEdges_newRep == afterStartIdx));     

rr = row1(ismember(row1,startEdgeIdx)==0);      % row1两个元素中不等于startIdx中的那个
rm = col1(ismember(row1,startEdgeIdx)==0);      % 1或者2，上面选取的索引对应的位置。

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

if s_b * s_t <= 0   % 11hit  不同号 
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
       endIdx = s(end);        % roww中最后一个满足小于length(patientEdge)/2的索引；
       addTris = [addTris;[patientEdge(roww(endIdx):length(patientEdge),:),repmat(rootEdges(endIdx,2),length(patientEdge)-roww(endIdx)+1,1)]];
       addTris = [addTris;[patientEdge(1:roww(end),:),repmat(rootEdges(endIdx,2),roww(end),1)]]; 
       addTris = [addTris;[patientEdge(roww(end):roww(1)-1,:),repmat(rootEdges(1,1),roww(1)-roww(end),1)]]; 
    end

else       %12hit
    s = find(roww < length(patientEdge)/2);      
    start = s(1);            % roww中第一个满足小于length(patientEdge)/2的索引；
   
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

writeOBJ('变形前的合并网格.obj', allVers, newTris);
OBJwriteVertices('centerPatient.obj', centerPatient);


%% 5  变形控制部分
index_crown_change = find(patientTooth.vertex(:,2)<(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+1 ...
                        & patientTooth.vertex(:,2)>=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-2);
patientTransform = patientTooth.vertex(index_crown_change,:);

%牙冠牙根分别形成网格
mergedToothVers = allVers;
mergeRegionIdx =  find(mergedToothVers(:,2)<(centerPatient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+1 ...
    &mergedToothVers(:,2)>=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-1);
mergeRegionVers = mergedToothVers(mergeRegionIdx,:);

% 病人牙冠部分多截一部分出来，用以确定标准牙根变形后点的位置，找到movedRootTooth2中离crownforchange中最近点的位置，并将这些点变为新的点，即牙冠点的位置
indices = 1:size(mergedToothVers,1);

% 合并网格中，不需要参与变形的顶点的索引，是一个行向量。
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
 

%% 6. 变形
[omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(mergedToothVers,1), newTris, exterior);

%[bi_L,bi_U,bi_P,bi_Q,bi_D,bi_S,bi_M] = biharm_factor_system(mergedToothVers,newTris, 'ext','voronoi', 'no_flatten',omega,N0,N1);
[A, bi_S] =  biharm_factor_system_modified(mergedToothVers,newTris,omega,N0,N1);

%找到牙根部分跟牙冠变形部分最近的点，将牙根上的点拉至牙冠处
finalVers = mergedToothVers;
for i = 1:length(mergeRegionVers)
   [minValue,startEdgeIdx] = mindis(patientTransform, mergeRegionVers(i,:),1);
   minvalue(i) = minValue;  row(i) = startEdgeIdx;
   finalVers(mergeRegionIdx(i),:) = patientTransform(startEdgeIdx,:);
end

 
BZ1 = zeros(size(mergedToothVers,1),3);
finalVers = biharm_solve_with_factor_modified(A, bi_S, finalVers, omega, N0, N1);
 
writeOBJ('最终网格.obj',finalVers, newTris)


%% 7. 纠正三角片方向
focusVers = finalVers(mergeRegionIdx, :);
OBJwriteVertices('focusVers.obj', focusVers);

focusCenter = mean(focusVers);      % 牙齿拼接处附近所有顶点的中点。
deltaHeight = 3;


disp('finished.');
