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




% 2.
firstVerIdx = newRootEdges(:,1);
 
secondVerIdx = newRootEdges(:,2);
 
middle = (allVers(firstVerIdx, :) + allVers(secondVerIdx, :))/2.0;      % ������Ե�ߵ��е�

pEdgeVers = allVers(patientEdgeVersIdx,:);
tempK = createns(pEdgeVers,'nsmethod','kdtree');
nearestIdx = knnsearch(tempK, middle, 'K', 1);                      % �ҳ��������ڱ�Ե���о���ÿ���е�����ĵ��������
addTris = [newRootEdges, patientEdgeVersIdx(nearestIdx)];       % ������Ե�����˵������������һ���µ�����Ƭ��

 
 

%% 3
if nearestIdx(1) < length(patientEdge)/2    % hit
    for j = 1:length(nearestIdx)-1
        
        if nearestIdx(j) ~= nearestIdx(j+1) 
            % nearestIdx(j)~nearestIdx(j+1)���������ڱ�Ե�ߺ�������Ե�ϵĵ�newRootEdges(j,2)���������Ƭ��
            temp1 = repmat(newRootEdges(j,2), nearestIdx(j+1)-nearestIdx(j), 1);
            tempa = nearestIdx(j);
            tempb = nearestIdx(j+1)-1;
            vecIdx = tempa:tempb;
            newTri = [patientEdge(vecIdx, :), temp1];
            addTris = [addTris; newTri];
            disp(vecIdx);
        end
    end
    
    lastIdx = nearestIdx(end);

    if lastIdx==length(patientEdge)
        newTri = [patientEdge(lastIdx,:),newRootEdges(1,1)];
        addTris = [addTris; newTri];
        
        if nearestIdx(1)~=1
            newTri = [patientEdge(1:nearestIdx(1),:),repmat(newRootEdges(1,1),nearestIdx(1),1)];
            addTris = [addTris; newTri];
        end

    elseif lastIdx > length(patientEdge)/2 && lastIdx<length(patientEdge)      
        newTri = [patientEdge(lastIdx:length(patientEdge),:),repmat(newRootEdges(1,1),length(patientEdge)-lastIdx,1)];
        addTris = [addTris; newTri];

        if nearestIdx(1)~=1
            addTris = [addTris;[patientEdge(1:nearestIdx(1),:),repmat(newRootEdges(1,1),nearestIdx(1),1)]];
        end

    elseif  lastIdx < length(patientEdge)/2 && nearestIdx(1)~=1      % hit
       plen = length(patientEdge);
       vecIdx = find(nearestIdx > plen/2);
       temp3 = vecIdx(end);
       temp4 = nearestIdx(temp3);
       temp1 = patientEdge(temp4: plen, :);
       temp2 = repmat(newRootEdges(temp3, 2), plen - temp4+1, 1);
       temp = [temp1,temp2];
       addTris = [addTris; temp];
       
       temp1 = patientEdge(1: lastIdx,:);
       temp2 = repmat(newRootEdges(temp3,2), lastIdx,1);
       temp = [temp1,temp2];
       addTris = [addTris; temp]; 
       
       
       temp1 = patientEdge(lastIdx: nearestIdx(1)-1,:);
       temp2 = repmat(newRootEdges(1,1),nearestIdx(1)-lastIdx,1);
       temp = [temp1, temp2];
       addTris = [addTris;temp]; 
    end

else  % pass
   vecIdx = find(nearestIdx < length(patientEdge)/2);
   for j = vecIdx(1):length(nearestIdx)-1
        if nearestIdx(j) ~= nearestIdx(j+1) 
            newTri = [patientEdge(nearestIdx(j):nearestIdx(j+1)-1,:), repmat(newRootEdges(j,2), nearestIdx(j+1)-nearestIdx(j),1)];
            addTris = [addTris; newTri];   
        end

   end 

   if vecIdx(1)>2
        for j = 1:vecIdx(1)-2
             if nearestIdx(j) ~= nearestIdx(j+1) 
                addTris = [addTris;[patientEdge(nearestIdx(j):nearestIdx(j+1)-1,:),repmat(newRootEdges(j,2),nearestIdx(j+1)-nearestIdx(j),1)]];   
             end
        end
   end

   if lastIdx < nearestIdx(1) 
      addTris = [addTris;[patientEdge(lastIdx:nearestIdx(1)-1,:),repmat(newRootEdges(end,2),nearestIdx(1)-lastIdx,1)]];
   end

   if nearestIdx(vecIdx(1)-1) <= length(patientEdge)
         addTris = [addTris;[patientEdge(nearestIdx(vecIdx(1)-1):length(patientEdge),:),repmat(newRootEdges(vecIdx(1)-1,2),length(patientEdge)-nearestIdx(vecIdx(1)-1)+1,1)]];
   end

   addTris = [addTris;[patientEdge(1:nearestIdx(vecIdx(1))-1,:),repmat(newRootEdges(vecIdx(1),1),nearestIdx(vecIdx(1))-1,1)]];
end


newTris = [allTris; addTris];

disp('finished.');
