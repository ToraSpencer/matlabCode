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
[~,r] = mindis(allVers(rootBdryEdges_newRep(:,1),:), allVers(bdryEdges_newRep(1,1),:), 1);


% 判断[edg_b(1,：),Edge_t(r,1)]de 方向与[Edge_t(r,：),edg_b(1,1)]的方向是否相同
% 为了防止个别特殊情况，找Edge_t(r,1)后面的第二个点
[r2,r3] = ind2sub(size(rootBdryEdges_newRep),find(rootBdryEdges_newRep == rootBdryEdges_newRep(r,2)));
rr = r2(ismember(r2,r)==0);
rm = r3(ismember(r2,r)==0);
tri11 = allVers(patientEdge(3,1),:) - allVers(patientEdge(1,1),:);
tri12 = allVers(rootBdryEdges_newRep(r,1),:) - allVers(patientEdge(3,1),:);
tri21 = allVers(rootBdryEdges_newRep(rr,3-rm),:) - allVers(rootBdryEdges_newRep(r,1),:);
tri22 = allVers(patientEdge(1,1),:) - allVers(rootBdryEdges_newRep(rr,3-rm),:);

s_b = sign(dot(allVers(patientEdge(1,1),:) -mean(rootEdgeVers),cross(tri11,tri12)));
s_t = sign(dot(allVers(rootBdryEdges_newRep(r,1),:) - mean(rootEdgeVers),cross(tri21,tri22)));


%% 对根部分的边界点排序，要求？？？
if s_b ~= s_t   % hit
    edg_t = sotr_edge(rootBdryEdges_newRep, r);
else
    E_t = [rootBdryEdges_newRep(:,2),rootBdryEdges_newRep(:,1)];
    edg_t = sotr_edge(E_t,r);
end


addTris = [];           % 补的三角片
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
