%% 测试4. 合并网格补三角片

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


%5.20更新


%1.所有的边界点按照一个方向排序
firstVer = bdryEdges_newRep(1,1);
[~,nearestIdx] = mindis(allVers(rootBdryEdges_newRep(:,1),:), allVers(firstVer,:), 1);


% 判断[edg_b(1,：),Edge_t(r,1)]de 方向与[Edge_t(r,：),edg_b(1,1)]的方向是否相同
% 为了防止个别特殊情况，找Edge_t(r,1)后面的第二个点
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


%% 1. 对根部分的边界点排序，要求？？？
if s_b ~= s_t   % hit
    newRootEdges = sotr_edge(rootBdryEdges_newRep, nearestIdx);
else   % pass
    E_t = [rootBdryEdges_newRep(:,2),rootBdryEdges_newRep(:,1)];
    newRootEdges = sotr_edge(E_t,nearestIdx);
end




% 2.
firstVerIdx = newRootEdges(:,1);
 
secondVerIdx = newRootEdges(:,2);
 
middle = (allVers(firstVerIdx, :) + allVers(secondVerIdx, :))/2.0;      % 牙根边缘边的中点

pEdgeVers = allVers(patientEdgeVersIdx,:);
tempK = createns(pEdgeVers,'nsmethod','kdtree');
nearestIdx = knnsearch(tempK, middle, 'K', 1);                      % 找出病人牙冠边缘点中距离每个中点最近的点的索引。
addTris = [newRootEdges, patientEdgeVersIdx(nearestIdx)];       % 牙根边缘边两端点和最近点组成了一个新的三角片。

 
 

%% 3
if nearestIdx(1) < length(patientEdge)/2    % hit
    for j = 1:length(nearestIdx)-1
        
        if nearestIdx(j) ~= nearestIdx(j+1) 
            % nearestIdx(j)~nearestIdx(j+1)的所有牙冠边缘边和牙根边缘上的点newRootEdges(j,2)组成新三角片。
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
