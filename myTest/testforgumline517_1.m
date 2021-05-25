clc
clear all
%%
functionname='testforgumline517_1.m'; 
functiondir=which(functionname);
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

%%  1.

% % 加载标准牙
load('dental_crown.mat');
load('dentalmodelwithroot0.1forlow.mat');
load('axisofdental_crowm');
load('trfromtoothtocrown.mat');
 

x = 11;             % 取11号牙来测试
for i = 1:28
    if  dentalwithtooth1(i).ID == x
        rootTooth = dentalwithtooth1(i);         % 有牙根的标准牙网格
    end
end
 
for j = 1:28
    if (upax(j).fid == x)
        axisStandard = upax(j);
    end
end

for k = 1:28
    if (dental_crown(k).fid == x)
        crown = dental_crown(k);            % 只有牙冠的标准牙网格？？？
    end
end


s = fix(x/10);  %取整
g = mod(x,10); %取余

fdi = textread('FDIUpper__.dxt');
toothIdx = find (fdi == x);
namestr1 = ['toothUpper_',num2str(toothIdx-1),'.','obj'];
patientTooth = Read_Obj(namestr1);                       % 病人的牙冠网格
namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];

gumline =  ReadObj(namestr2);
axis = ReadObj('AXISUpper_.obj');
axisPatient = axis( (3*(toothIdx-1) + 1) : 3*toothIdx, :);       % 病人牙齿的三个牙轴方向向量
crownTooth = crown.model;          % 标准牙冠



%       1.2 确定病人牙冠的切割部分
centerPatient = mean(patientTooth.vertex);
y = max(gumline(:,2));
p_patient =[centerPatient(1),y(1),centerPatient(3)];
patientYdir = axisPatient(2,:);
 

cutPatientIdx = find(patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-2);
cutPatientCrownVers = patientTooth.vertex(cutPatientIdx,:);



%       1.3 确定病人牙冠边界
a1 = ismember(patientTooth.face(:,1),cutPatientIdx);
aa1 = find(a1==1);
a2 = ismember(patientTooth.face(:,2),cutPatientIdx);
aa2 = find(a2==1);
a3 = ismember(patientTooth.face(:,3),cutPatientIdx);
aa3 = find(a3==1);
[c1,~,~] = intersect(aa1,aa2);c2 = intersect(c1,aa3);


%找边界
[raw_edges_list] = query_edges_list(patientTooth.face(c2,:),'sorted');
[~,~,iu] = unique(sort(raw_edges_list,2),'rows');
[i3,i4] = histc(iu,unique(iu));
lone_edges_idx_vect = i3(i4) == 1;
lone_edges_list = unique(raw_edges_list(lone_edges_idx_vect,:),'rows');

% 边界点：
edge_p_b_local_list = unique([lone_edges_list(:,1);lone_edges_list(:,2)]);

[C_b,edge_p_b_local,~] = intersect(cutPatientIdx,edge_p_b_local_list);
edge_p_b = cutPatientCrownVers(edge_p_b_local,:);
Edge_b = zeros(size(lone_edges_list));
for i = 1:length(edge_p_b_local)
    nu = find(lone_edges_list == edge_p_b_local_list(i));
    Edge_b(nu) = edge_p_b_local(i);
end



patientEdge = sotr_edge(Edge_b,1);


ff_b = patientTooth.face(c2,:);
patientEdgeTris = ones (size(ff_b));
for j = 1:length(cutPatientCrownVers)
    nu = find( ff_b == cutPatientIdx(j));
    patientEdgeTris(sub2ind(size(patientEdgeTris), nu)) = j;
end



%        1.4 确定带根标准牙切割参数
centerCrown = mean(crownTooth.vertex);
rootYdir = axisStandard.y;
index_crown_biaozhun = find(crownTooth.vertex(:,2)<=(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(:,1) ...
            - centerCrown(1))+rootYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/rootYdir(2)));
v_crown_biaozhun = abs(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(index_crown_biaozhun,1)-centerCrown(1))...
             +rootYdir(3)*(crownTooth.vertex(index_crown_biaozhun,3)-centerCrown(3)))/rootYdir(2));
in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
point_crown =  crownTooth.vertex(index_crown_biaozhun(in_crown),:);

%标准牙根
p_tooth_biaozhun = mean(rootTooth.vertices);
index_tooth_biaozhun = find(rootTooth.vertices(:,2)<=(p_tooth_biaozhun(2) - (rootYdir(1)*(rootTooth.vertices(:,1) ...
                 - p_tooth_biaozhun(1))+rootYdir(3)*(rootTooth.vertices(:,3) - p_tooth_biaozhun(3)))/rootYdir(2)));
vz_tooth = abs(p_tooth_biaozhun(2) - (rootYdir(1)*(rootTooth.vertices(index_tooth_biaozhun,1)-p_tooth_biaozhun(1))...
           +rootYdir(3)*(rootTooth.vertices(index_tooth_biaozhun,3)-p_tooth_biaozhun(3)))/rootYdir(2));
in_tooth = find(vz_tooth == max(vz_tooth));
point_tooth = rootTooth.vertices(index_tooth_biaozhun(in_tooth),:);
centcrownintooth = centerCrown+point_tooth-point_crown;%切后标准根的边缘中心



%%
% 2. 对齐——利用标准牙和病人牙齿的中心点以及三轴进行对齐
R = inv([axisStandard.x;axisStandard.y;axisStandard.z])*axisPatient;
C=centcrownintooth*R;
tooth_T = (rootTooth.vertices) *R + repmat((centerPatient - C),length(rootTooth.vertices),1);

%       2.1 旋转平移带根标准牙
index_1 = find(patientTooth.vertex(:,2)>=(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))...
        & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+0.5);
point1 = patientTooth.vertex(index_1,:);
index_2 = find(tooth_T(:,2)>=(centerPatient(2) - (patientYdir(1)*(tooth_T(:,1) - p_patient(1))+patientYdir(3)*(tooth_T(:,3) - p_patient(3)))/patientYdir(2))...
        & tooth_T(:,2)<=(p_patient(2) - (patientYdir(1)*(tooth_T(:,1) - p_patient(1))+patientYdir(3)*(tooth_T(:,3) - p_patient(3)))/patientYdir(2))+0.5);
point2 = tooth_T(index_2,:);

for i = 1:length(point2)
   [minValue,r]=mindis(point1,point2(i,:),1);
   minvalue(i) = minValue;  row(i) = r;
end

point11 = point1(row,:);
[R,t,BRt,e,~,~] = icp(point11,point2);
tooth_T = bsxfun(@plus,tooth_T*R,t);

save('tooth_T.mat', 'tooth_T');
OBJwriteVertices('第二次旋转平移后的带根标准牙.obj', tooth_T);



%%
% 3. 切割——确定带根标准牙的切割部分。

%       3.1 

rootCutIdx =  find(tooth_T(:,2)>(p_patient(2) - ( patientYdir(1)*(tooth_T(:,1) ...
                 - p_patient(1))+ patientYdir(3)*(tooth_T(:,3) - p_patient(3)))/ patientYdir(2))-2);%切出来的牙根的索引值%注意！按照病人牙冠方向切
rootCutVers = tooth_T(rootCutIdx ,:);


OBJwriteVertices('合并网格的牙根部分点集.obj', rootCutVers);

%% 4. 处理三角片信息

a_t1 = ismember(rootTooth.faces(:,1),rootCutIdx);
aa_t1 = find(a_t1==1);
a_t2 = ismember(rootTooth.faces(:,2),rootCutIdx);
aa_t2 = find(a_t2==1);
a_t3 = ismember(rootTooth.faces(:,3),rootCutIdx);
aa_t3 = find(a_t3==1);
[c_t1,~,~] = intersect(aa_t1,aa_t2);c_t2 = intersect(c_t1,aa_t3);


[raw_edges_list_t] = query_edges_list(rootTooth.faces(c_t2,:),'sorted');
[~,~,iu_t] = unique(sort(raw_edges_list_t,2),'rows');
[i3,i4] = histc(iu_t,unique(iu_t));
lone_edges_idx_vect_t = i3(i4) == 1;
lone_edges_list_t = unique(raw_edges_list_t(lone_edges_idx_vect_t,:),'rows');


ff_t = rootTooth.faces(c_t2,:);
f_t0 = ones (size(ff_t));
for j = 1:length(rootCutVers)
    nu = find( ff_t == rootCutIdx(j));
    f_t0(sub2ind(size(f_t0), nu)) = j+length(cutPatientCrownVers);
end
f_end = [patientEdgeTris;f_t0];
v_end = [cutPatientCrownVers;rootCutVers];

%边界点：
edge_p_t_local_list = unique([lone_edges_list_t(:,1);lone_edges_list_t(:,2)]);

%边界点在牙根点中的位置
[C_t,edge_p_t_local,~] = intersect(rootCutIdx,edge_p_t_local_list);

edge_p_t = rootCutVers(edge_p_t_local,:);
edge_p_t_local = edge_p_t_local+length(cutPatientCrownVers);
Edge_t = zeros(size(lone_edges_list_t));

for i = 1:length(edge_p_t_local)
    nu = find(lone_edges_list_t == edge_p_t_local_list(i));
    Edge_t(nu) = edge_p_t_local(i);
end


%5.20更新
%1.所有的边界点按照一个方向排序
[~,r]=mindis(v_end(Edge_t(:,1),:),v_end(Edge_b(1,1),:),1);


%判断[edg_b(1,：),Edge_t(r,1)]de 方向与[Edge_t(r,：),edg_b(1,1)]的方向是否相同
%为了防止个别特殊情况，找Edge_t(r,1)后面的第二个点
[r2,r3] = ind2sub(size(Edge_t),find(Edge_t == Edge_t(r,2)));
rr = r2(ismember(r2,r)==0);rm = r3(ismember(r2,r)==0);
tri11 = v_end(patientEdge(3,1),:) - v_end(patientEdge(1,1),:);
tri12 = v_end(Edge_t(r,1),:) - v_end(patientEdge(3,1),:);
tri21 = v_end(Edge_t(rr,3-rm),:) - v_end(Edge_t(r,1),:);
tri22 = v_end(patientEdge(1,1),:) - v_end(Edge_t(rr,3-rm),:);
s_b = sign(dot(v_end(patientEdge(1,1),:) -mean(edge_p_b),cross(tri11,tri12)));
s_t = sign(dot(v_end(Edge_t(r,1),:) - mean(edge_p_t),cross(tri21,tri22)));


if s_b ~= s_t
   edg_t = sotr_edge(Edge_t,r);
else
    E_t = [Edge_t(:,2),Edge_t(:,1)];
    edg_t = sotr_edge(E_t,r);
end 

%边界点融合,以稀疏（牙根）的边线为基准，找密集（牙冠）中离两个点距离最近的点

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



%% 5  变形控制部分
index_crown_change = find(patientTooth.vertex(:,2)>(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-1 ...
                        & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+2);
patientTransform = patientTooth.vertex(index_crown_change,:);


%牙冠牙根分别形成网格
mergedToothVers = v_end;
index_tooth_change =  find(mergedToothVers(:,2)>(centerPatient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-1 ...
    &mergedToothVers(:,2)<=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+1);
mergeRegionVers = mergedToothVers(index_tooth_change,:);



%    合并网格中不需要变形的部分
indices = 1:size(mergedToothVers,1);
exterior = indices(mergedToothVers(:,2)<(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-0.6...
    |mergedToothVers(:,2)>=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+3);

 
%%
% 6. 变形
[omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(mergedToothVers,1), newTris, exterior);
[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(mergedToothVers,newTris, 'ext','voronoi', 'no_flatten',omega,N0,N1);


%找到牙根部分跟牙冠变形部分最近的点，将牙根上的点拉至牙冠处
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
 

% axis image
writeOBJ('最终网格.obj',finalVers, newTris)


disp('finished.');
