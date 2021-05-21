clc
close all;
clear all;
% functionname='test1.m'; 
% functiondir=which(functionname);
% functiondir=functiondir(1:end-length(functionname));
% addpath([functiondir 'mesh process'])
addpath('mesh process');
addpath('mixFE');

disp('start')


%% 1. 加载标准牙
load('dental_crown.mat');                       % 标准牙，只有牙冠没有牙根
load('dentalmodelwithroot0.1forlow.mat');       % 标准牙，有牙根
load('axisofdental_crowm');                     %            
 

%       1.1 相同的FDI号拿出来
x = 11;             % 取11号牙来测试
for i = 1:28
    if  dentalwithtooth1(i).ID == x
        rootTooth = dentalwithtooth1(i);         % 有牙根的标准牙网格
    end
end
 
for j = 1:28
    if (upax(j).fid == x)
        crown_ax = upax(j);
    end
end

for k = 1:28
    if (dental_crown(k).fid == x)
        crown = dental_crown(k);            % 只有牙冠的标准牙网格？？？
    end
end


s = fix(x/10);  %取整
g = mod(x,10); %取余


%% for debug
crownTooth = crown.model;
writeOBJ('标准牙冠.obj', crownTooth.vertex, crownTooth.face);
writeOBJ('带牙根标准牙.obj', rootTooth.vertices, rootTooth.faces);        % ！！！带牙根的标准牙网格三角片是反的，是否有问题？
crown_center = mean(crownTooth.vertex);
crown_ax_xline = getDirLine(crown_ax.x, crown_center);
crown_ax_yline = getDirLine(crown_ax.y, crown_center);
crown_ax_zline = getDirLine(crown_ax.z, crown_center);
OBJwriteVertices('标准牙x轴.obj', crown_ax_xline);
OBJwriteVertices('标准牙y轴.obj', crown_ax_yline);
OBJwriteVertices('标准牙z轴.obj', crown_ax_zline);
crown_ax_dirs = [crown_ax.x; crown_ax.y; crown_ax.z];
OBJwriteVertices('标准牙牙轴方向向量.obj', crown_ax_dirs); 
    fdi = textread('FDIUpper__.dxt');
    toothIdx = find (fdi == x);
    
    
namestr1 = ['toothUpper_',num2str(toothIdx-1),'.','obj'];
patientTooth = Read_Obj(namestr1);                       % 病人的牙冠网格
namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];
gumline =  ReadObj(namestr2);                       % 牙龈线点云
axis = ReadObj('AXISUpper_.obj');
axi = axis( (3*(toothIdx-1) + 1) : 3*toothIdx, :);       % 当前牙齿的三个牙轴方向向量
OBJwriteVertices('病人牙轴方向向量.obj', axi);
OBJwriteVertices('病人牙龈线.obj', gumline);


%       1.2 确定病人牙冠的切割部分
centerPatient = mean(patientTooth.vertex);     % 病人牙冠网格中心
ymax = max(gumline(:,2));             % 牙龈线点云中y坐标最大值
p_patient =[centerPatient(1), ymax, centerPatient(3)];
patientYdir = axi(2,:);                    % 当前病人牙齿的y轴方向
index_bingren = find(patientTooth.vertex(:,2) <= (p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) ...
        - p_patient(1)) + patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-2);    % 超参数为2
cutPatientVers = patientTooth.vertex(index_bingren,:);      % 牙齿自身坐标系下，只保留y方向上高于p_bingren点2mm的点。

disp(size(index_bingren));      % for debug


% FOR DEBUG
OBJwriteVertices('p_bingren.obj', p_patient);
moved_p_bingren = p_patient + 2 * patientYdir;
writeOBJ('病人牙冠网格.obj', patientTooth.vertex, patientTooth.face);
OBJwriteVertices('moved_p_bingren.obj', moved_p_bingren);
OBJwriteVertices('合并网格的病人牙冠部分点集.obj', cutPatientVers);
bingren_xline = getDirLine(axi(1,:), centerPatient);
bingren_yline = getDirLine(axi(2,:), centerPatient);
bingren_zline = getDirLine(axi(3,:), centerPatient);
OBJwriteVertices('病人牙x轴.obj', bingren_xline);
OBJwriteVertices('病人牙y轴.obj', bingren_yline);
OBJwriteVertices('病人牙z轴.obj', bingren_zline);
OBJwriteVertices('病人牙重心.obj', centerPatient);

centerCrown = mean(crownTooth.vertex);     % 标准牙冠重心
rootYdir = crown_ax.y;
index_crown = find(crownTooth.vertex(:,2)<=(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(:,1) ...
            - centerCrown(1))+rootYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/rootYdir(2)));
v_crown_biaozhun = abs(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(index_crown,1)-centerCrown(1))...
             +rootYdir(3)*(crownTooth.vertex(index_crown,3)-centerCrown(3)))/rootYdir(2));
in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
point_crown =  crownTooth.vertex(index_crown(in_crown),:);

%  for debug
OBJwriteVertices('point_crown.obj', point_crown);

%        1.3 确定带根标准牙切割部分
centerRoot = mean(rootTooth.vertices);
index_root_higher = find(rootTooth.vertices(:,2)<=(centerRoot(2) - (rootYdir(1)*(rootTooth.vertices(:,1) ...
                 - centerRoot(1))+rootYdir(3)*(rootTooth.vertices(:,3) ...
                 - centerRoot(3)))/rootYdir(2)));       % 牙齿自身坐标系下，y轴高度高于重心的点。
rootZvalue = abs(centerRoot(2) - (rootYdir(1)*(rootTooth.vertices(index_root_higher,1)-centerRoot(1))...
             +rootYdir(3)*(rootTooth.vertices(index_root_higher,3)-centerRoot(3)))/rootYdir(2));
in_tooth = find(rootZvalue == max(rootZvalue));
point_root = rootTooth.vertices(index_root_higher(in_tooth),:);
centCrownInTooth = centerCrown + point_root - point_crown;  % 切后标准根的边缘中心

% for debug
OBJwriteVertices('centerRoot.obj', centerRoot);  
OBJwriteVertices('point_root.obj', point_root);  
OBJwriteVertices('vz_tooth.obj', rootZvalue);  
OBJwriteVertices('标准牙中y坐标高于重心的点.obj', rootTooth.vertices(index_root_higher, :));
OBJwriteVertices('切后标准根的边缘中心.obj', centCrownInTooth);  




%%
% 2. 对齐——利用标准牙和病人牙齿的中心点以及三轴进行对齐
R = inv([crown_ax.x; crown_ax.y; crown_ax.z]) * axi;   % 旋转     

offset = centerPatient - centCrownInTooth * R;

tooth_T = (rootTooth.vertices) *R + repmat(offset,length(rootTooth.vertices),1);
new_crown_ax_xline = crown_ax_xline * R + repmat(offset,length(crown_ax_xline),1);
new_crown_ax_yline = crown_ax_yline * R + repmat(offset,length(crown_ax_yline),1);
new_crown_ax_zline = crown_ax_zline * R + repmat(offset,length(crown_ax_zline),1);
OBJwriteVertices('旋转平移后的标准牙x轴.obj', new_crown_ax_xline);
OBJwriteVertices('旋转平移后的标准牙y轴.obj', new_crown_ax_yline);
OBJwriteVertices('旋转平移后的标准牙z轴.obj', new_crown_ax_zline);
OBJwriteVertices('旋转平移后的标准牙重心.obj', mean(tooth_T));
OBJwriteVertices('centerRoot_new.obj', centerRoot + offset);

%       2.1 旋转平移带根标准牙
index_1 = find(patientTooth.vertex(:,2)>=(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))...
         & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - ...
         p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+0.5);
patientMergeVers = patientTooth.vertex(index_1,:);

index_2 = find(tooth_T(:,2)>=(centerPatient(2) - (patientYdir(1)*(tooth_T(:,1) - p_patient(1))+patientYdir(3)*(tooth_T(:,3) - p_patient(3)))/patientYdir(2))...
        & tooth_T(:,2)<=(p_patient(2) - (patientYdir(1)*(tooth_T(:,1) - p_patient(1))+patientYdir(3)*(tooth_T(:,3) - p_patient(3)))/patientYdir(2))+0.5);
rootMergeVers = tooth_T(index_2,:);

% for debug
 OBJwriteVertices('第一次旋转平移后的带根标准牙.obj', tooth_T);
 OBJwriteVertices('2.1病人牙冠上融合区域的点.obj', patientMergeVers);
 OBJwriteVertices('2.1带根标准牙融合区域的点.obj', rootMergeVers);

 
%       2.2 icp算法进一步对齐 
for i = 1:length(rootMergeVers)
   [minValue,r] = mindis(patientMergeVers, rootMergeVers(i,:),1);
   minvalue(i) = minValue;  row(i) = r;
end
point11 = patientMergeVers(row,:);
[R,t,BRt,e,~,~] = icp(point11, rootMergeVers);          % icp算法
tooth_T = bsxfun(@plus,tooth_T*R, t);   
OBJwriteVertices('第二次旋转平移后的带根标准牙.obj', tooth_T);
       


%%
% 3. 切割——对齐后按照病人牙冠的切割位置进行切割

%       3.1 确定带根标准牙的切割部分。
index_root_cut =  find(tooth_T(:,2)>(p_patient(2) - ( patientYdir(1)*(tooth_T(:,1) ...
                 - p_patient(1))+ patientYdir(3)*(tooth_T(:,3) - p_patient(3)))/ patientYdir(2))-2);  %切出来的牙根的索引值%注意！按照病人牙冠方向切
rootCutVers = tooth_T(index_root_cut ,:);


%       3.2 找出牙冠变形时需要参考的原病人牙齿网格上的顶点。
index_patientTransform = find(patientTooth.vertex(:,2)>(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-1 ...
                        & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+2);
patientTransform = patientTooth.vertex(index_patientTransform,:);

% for debug
OBJwriteVertices('合并网格的牙根部分点集.obj', rootCutVers);
OBJwriteVertices('合并区域参考的病人牙冠顶点.obj', patientTransform);



%%
% 4. 合并网格
%       4.1 牙冠和牙根合并成一个网格网格 
mergedToothVers = [cutPatientVers; rootCutVers];
[mergedToothFace]=MyCrustOpen(mergedToothVers);
% for debug
writeOBJ('合并网格.Obj', mergedToothVers, mergedToothFace);             % ！！！ 三角片朝向混乱



%%
% 5. 补洞
mergedToothBdr = select_holes_and_boundary(mergedToothVers, mergedToothFace);
newTris = fill_mesh_holes(mergedToothVers, mergedToothFace, mergedToothBdr,'closed',200);
newTris = double(newTris);
mergeRegionIdx =  find(mergedToothVers(:,2)>(centerPatient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-1 ...
    &mergedToothVers(:,2)<=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+1);


mergeRegionVers = mergedToothVers(mergeRegionIdx,:);



newTris = reduceWrongTris(newTris);     % 去除面信息中错误的三角片。


writeOBJ('补洞后的网格.obj', mergedToothVers, newTris);
OBJwriteVertices('mergeRegionVers.obj', mergeRegionVers);
vecWriteUnsigned('mergeRegionIdx.obj', mergeRegionIdx);



%    病人牙冠部分多截一部分出来，用以确定标准牙根变形后点的位置，找到tooth_T中离patientTransform中最近点的位置，并将这些点变为新的点，即牙冠点的位置
indices = 1:size(mergedToothVers,1);

%    合并网格中不需要变形的部分
exterior = indices(mergedToothVers(:,2)<(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-0.6...
    |mergedToothVers(:,2)>=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+3);
 


%%
% 6. 变形

% FOR DEBUG: 打印输入参数：
OBJwriteVertices('合并区域的顶点.obj', mergeRegionVers);
vecWriteUnsigned('非合并区域的顶点索引向量.obj', exterior);
save('非合并区域的顶点索引向量.mat', 'exterior');
OBJwriteVertices('合并区域参考的病人牙冠顶点.obj', patientTransform);
writeOBJ('合并网格', mergedToothVers, newTris);



%       6.1 计算变形所需的参数
[omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(mergedToothVers,1), newTris, exterior);
save('mergedToothVers.mat', 'mergedToothVers');
save('newTris.mat', 'newTris');
save('omega.mat', 'omega');
save('N0.mat','N0');
save('N1.mat','N1');
save('mergeRegionVers.mat','mergeRegionVers');
save('patientTransform.mat', 'patientTransform');
save('mergeRegionIdx.mat', 'mergeRegionIdx');
%[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(mergedToothVers, newTris, 'ext', 'voronoi', 'no_flatten',omega, N0, N1);

[bi_L,bi_U,bi_P,bi_Q,bi_D,bi_S,bi_M] = biharm_factor_system_modified(mergedToothVers, newTris,omega, N0, N1);



%       6.2 执行变形——找到牙根部分跟牙冠变形部分最近的点，将牙根上的点拉至牙冠处
finalVers = mergedToothVers;
for i = 1:length(mergeRegionVers)
   [minValue,r]=mindis(patientTransform, mergeRegionVers(i,:), 1);
   minvalue(i) = minValue;  row(i) = r;
   finalVers(mergeRegionIdx(i),:) = patientTransform(r,:);
end

save('bi_L.mat', 'bi_L');
save('bi_U.mat', 'bi_U');
save('bi_P.mat', 'bi_P');
save('bi_Q.mat', 'bi_Q');
save('bi_D.mat', 'bi_D');
save('bi_S.mat', 'bi_S');
save('bi_M.mat', 'bi_M');

% finalVers = biharm_solve_with_factor( ...
%     bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
%     newTris, finalVers, omega, N0, N1, 'ext', 'no_flatten', BZ1, mergedToothVers);

finalVers = biharm_solve_with_factor_modified( ...
    bi_L, bi_U, bi_P, bi_Q, bi_D, bi_S, ...
     finalVers, omega, N0, N1);


writeOBJ('最终结果.obj', finalVers, newTris);


disp('finished');

 
