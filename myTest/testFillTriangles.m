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
addTris = [];           % 补的三角片

tempMiddle = [];
for j = 1:length(newRootEdges)
    p1 = allVers(newRootEdges(j,1),:);      % 牙根边缘第j条边第一个点
    p2 = allVers(newRootEdges(j,2),:);      % 

    middle = (p1 + p2)/2;    
    tempMiddle = [tempMiddle; middle];
    [~,nearestIdx] = mindis(allVers(patientEdgeVersIdx,:), middle, 1);
    minDisIdx(j) = nearestIdx;              % 牙根边缘中，距离第j条边中点最近的点的索引。
    newTri = [newRootEdges(j,:), patientEdgeVersIdx(nearestIdx)];
    addTris = [addTris; newTri];    % 找出的最近点所在的边的两端点，和当前的边的前一个端点构成了一个补的三角片。
end


%% 3
if minDisIdx(1) < length(patientEdge)/2    % hit
    for j = 1:length(minDisIdx)-1
        
        if minDisIdx(j) ~= minDisIdx(j+1) 
            % minDisIdx(j)~minDisIdx(j+1)的所有牙冠边缘边和牙根边缘上的点newRootEdges(j,2)组成新三角片。
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
