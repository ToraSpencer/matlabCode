%% ����4. �ϲ���������Ƭ

clc;
clear all;
close all;
%%
load('rootBdryEdges_newRep.mat');
load('patientEdge.mat');
load('allVers.mat');
load('allTris.mat');
load('rootEdgeVers.mat');
load('bdryEdges_newRep.mat');


%5.20����


%1.���еı߽�㰴��һ����������
firstVer = bdryEdges_newRep(1,1);
[~,nearestIdx] = mindis(allVers(rootBdryEdges_newRep(:,1),:), allVers(firstVer,:), 1);


% �ж�[edg_b(1,��),Edge_t(r,1)]de ������[Edge_t(r,��),edg_b(1,1)]�ķ����Ƿ���ͬ
% Ϊ�˷�ֹ���������������Edge_t(r,1)����ĵڶ�����
[r2,r3] = ind2sub(size(rootBdryEdges_newRep),find(rootBdryEdges_newRep == rootBdryEdges_newRep(nearestIdx,2)));
rr = r2(ismember(r2,nearestIdx)==0);
rm = r3(ismember(r2,nearestIdx)==0);
tri11 = allVers(patientEdge(3,1),:) - allVers(patientEdge(1,1),:);
tri12 = allVers(rootBdryEdges_newRep(nearestIdx,1),:) - allVers(patientEdge(3,1),:);
tri21 = allVers(rootBdryEdges_newRep(rr,3-rm),:) - allVers(rootBdryEdges_newRep(nearestIdx,1),:);
tri22 = allVers(patientEdge(1,1),:) - allVers(rootBdryEdges_newRep(rr,3-rm),:);

s_b = sign(dot(allVers(patientEdge(1,1),:) -mean(rootEdgeVers),cross(tri11,tri12)));
s_t = sign(dot(allVers(rootBdryEdges_newRep(nearestIdx,1),:) - mean(rootEdgeVers),cross(tri21,tri22)));

patientEdgeVersIdx = patientEdge(:,1);


%% 1. �Ը����ֵı߽������Ҫ�󣿣���
if s_b ~= s_t   % hit
    newRootEdges = sotr_edge(rootBdryEdges_newRep, nearestIdx);
else   % pass
    E_t = [rootBdryEdges_newRep(:,2),rootBdryEdges_newRep(:,1)];
    newRootEdges = sotr_edge(E_t,nearestIdx);
end




%% 2.
 
% firstVerIdx = newRootEdges(:,1);
%  
% secondVerIdx = newRootEdges(:,2);
%  
% middle = (allVers(firstVerIdx, :) + allVers(secondVerIdx, :))/2.0;
% 
% tempK = createns(allVers(patientEdgeVersIdx,:),'nsmethod','kdtree');
% nearestIdx = knnsearch(tempK, middle, 'K', 1); 
% addTris = [newRootEdges, patientEdgeVersIdx(nearestIdx)];

 

%% 2.
addTris = [];           % ��������Ƭ

tempMiddle = [];
for j = 1:length(newRootEdges)
    p1 = allVers(newRootEdges(j,1),:);      % ������Ե��j���ߵ�һ����
    p2 = allVers(newRootEdges(j,2),:);      % 

    middle = (p1 + p2)/2;    
    tempMiddle = [tempMiddle; middle];
    [~,nearestIdx] = mindis(allVers(patientEdgeVersIdx,:), middle, 1);
    minDisIdx(j) = nearestIdx;              % ������Ե�У������j�����е�����ĵ��������
    newTri = [newRootEdges(j,:), patientEdgeVersIdx(nearestIdx)];
    addTris = [addTris; newTri];    % �ҳ�����������ڵıߵ����˵㣬�͵�ǰ�ıߵ�ǰһ���˵㹹����һ����������Ƭ��
end


%% 3
if minDisIdx(1) < length(patientEdge)/2    % hit
    for j = 1:length(minDisIdx)-1
        
        if minDisIdx(j) ~= minDisIdx(j+1) 
            % minDisIdx(j)~minDisIdx(j+1)���������ڱ�Ե�ߺ�������Ե�ϵĵ�newRootEdges(j,2)���������Ƭ��
            temp = repmat(newRootEdges(j,2), minDisIdx(j+1)-minDisIdx(j), 1);
            newTri = [patientEdge(minDisIdx(j): minDisIdx(j+1)-1, :), temp];
            addTris = [addTris; newTri];
        end
    end

    if minDisIdx(end)==length(patientEdge)
        newTri = [patientEdge(minDisIdx(end),:),newRootEdges(1,1)];
        addTris = [addTris; newTri];
        
        if minDisIdx(1)~=1
            newTri = [patientEdge(1:minDisIdx(1),:),repmat(newRootEdges(1,1),minDisIdx(1),1)];
            addTris = [addTris; newTri];
        end

    elseif minDisIdx(end) > length(patientEdge)/2 && minDisIdx(end)<length(patientEdge)
        newTri = [patientEdge(minDisIdx(end):length(patientEdge),:),repmat(newRootEdges(1,1),length(patientEdge)-minDisIdx(end),1)];
        addTris = [addTris; newTri];

        if minDisIdx(1)~=1
            addTris = [addTris;[patientEdge(1:minDisIdx(1),:),repmat(newRootEdges(1,1),minDisIdx(1),1)]];
        end

    elseif  minDisIdx(end) < length(patientEdge)/2 && minDisIdx(1)~=1

       s = find(minDisIdx > length(patientEdge)/2);
       addTris = [addTris;[patientEdge(minDisIdx(s(end)):length(patientEdge),:),repmat(newRootEdges(s(end),2),length(patientEdge)-minDisIdx(s(end))+1,1)]];
       addTris = [addTris;[patientEdge(1:minDisIdx(end),:),repmat(newRootEdges(s(end),2),minDisIdx(end),1)]]; 
       addTris = [addTris;[patientEdge(minDisIdx(end):minDisIdx(1)-1,:),repmat(newRootEdges(1,1),minDisIdx(1)-minDisIdx(end),1)]]; 
    end

else  % pass
   s = find(minDisIdx < length(patientEdge)/2);
   for j = s(1):length(minDisIdx)-1
        if minDisIdx(j) ~= minDisIdx(j+1) 
            newTri = [patientEdge(minDisIdx(j):minDisIdx(j+1)-1,:), repmat(newRootEdges(j,2), minDisIdx(j+1)-minDisIdx(j),1)];
            addTris = [addTris; newTri];   
        end

   end 

   if s(1)>2
        for j = 1:s(1)-2
             if minDisIdx(j) ~= minDisIdx(j+1) 
                addTris = [addTris;[patientEdge(minDisIdx(j):minDisIdx(j+1)-1,:),repmat(newRootEdges(j,2),minDisIdx(j+1)-minDisIdx(j),1)]];   
             end
        end
   end

   if minDisIdx(end) < minDisIdx(1) 
      addTris = [addTris;[patientEdge(minDisIdx(end):minDisIdx(1)-1,:),repmat(newRootEdges(end,2),minDisIdx(1)-minDisIdx(end),1)]];
   end

   if minDisIdx(s(1)-1) <= length(patientEdge)
         addTris = [addTris;[patientEdge(minDisIdx(s(1)-1):length(patientEdge),:),repmat(newRootEdges(s(1)-1,2),length(patientEdge)-minDisIdx(s(1)-1)+1,1)]];
   end

   addTris = [addTris;[patientEdge(1:minDisIdx(s(1))-1,:),repmat(newRootEdges(s(1),1),minDisIdx(s(1))-1,1)]];
end



newTris = [allTris; addTris];

disp('finished.');
