clc
close all;
clear all;
% functionname='test1.m'; 
% functiondir=which(functionname);
% functiondir=functiondir(1:end-length(functionname));
% addpath([functiondir 'mesh process'])
addpath(['mesh process'])

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


if s ==1 || s ==2               % 上颌牙
    fdi = textread('FDIUpper__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
        disp( '该病例没有此牙齿 ');
        return;
       
    else
        namestr1 = ['toothUpper_',num2str(toothIdx-1),'.','obj'];
        patientTooth = Read_Obj(namestr1);                       % 病人的牙冠网格
        namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);                       % 牙龈线点云
        axis = ReadObj('AXISUpper_.obj');
        axi = axis( (3*(toothIdx-1) + 1) : 3*toothIdx, :);       % 当前牙齿的三个牙轴方向向量
        OBJwriteVertices('病人牙轴方向向量.obj', axi);
        OBJwriteVertices('病人牙龈线.obj', gumline);
        
        
        % 1.2 确定病人牙冠的切割部分
        centerPatient = mean(patientTooth.vertex);     % 病人牙冠网格中心
        ymax = max(gumline(:,2));             % 牙龈线点云中y坐标最大值
        p_bingren =[centerPatient(1), ymax, centerPatient(3)];
        toothYdir = axi(2,:);                    % 当前病人牙齿的y轴方向
        index_bingren = find(patientTooth.vertex(:,2) <= (p_bingren(2) - (toothYdir(1)*(patientTooth.vertex(:,1) ...
                - p_bingren(1)) + toothYdir(3)*(patientTooth.vertex(:,3) - p_bingren(3)))/toothYdir(2))-2);    % 超参数为2
        cutPatientVers = patientTooth.vertex(index_bingren,:);      % 牙齿自身坐标系下，只保留y方向上高于p_bingren点2mm的点。
        
        disp(size(index_bingren));      % for debug
        
        
       % FOR DEBUG
        OBJwriteVertices('p_bingren.obj', p_bingren);
        moved_p_bingren = p_bingren + 2 * toothYdir;
        writeOBJ('病人牙冠网格.obj', patientTooth.vertex, patientTooth.face);
        OBJwriteVertices('1.2moved_p_bingren.obj', moved_p_bingren);
        OBJwriteVertices('1.2除去底部的病人牙冠点集.obj', cutPatientVers);
        bingren_xline = getDirLine(axi(1,:), centerPatient);
        bingren_yline = getDirLine(axi(2,:), centerPatient);
        bingren_zline = getDirLine(axi(3,:), centerPatient);
        OBJwriteVertices('病人牙x轴.obj', bingren_xline);
        OBJwriteVertices('病人牙y轴.obj', bingren_yline);
        OBJwriteVertices('病人牙z轴.obj', bingren_zline);
        
        centerCrown = mean(crownTooth.vertex);     % 标准牙冠重心
        crownYdir = crown_ax.y;
        index_crown = find(crownTooth.vertex(:,2)<=(centerCrown(2) - (crownYdir(1)*(crownTooth.vertex(:,1) ...
                    - centerCrown(1))+crownYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/crownYdir(2)));
        v_crown_biaozhun = abs(centerCrown(2) - (crownYdir(1)*(crownTooth.vertex(index_crown,1)-centerCrown(1))...
                     +crownYdir(3)*(crownTooth.vertex(index_crown,3)-centerCrown(3)))/crownYdir(2));
        in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
        point_crown =  crownTooth.vertex(index_crown(in_crown),:);
        
       %  for debug
        OBJwriteVertices('point_crown.obj', point_crown);
        
       % 1.3 确定带根标准牙切割部分
        centerRoot = mean(rootTooth.vertices);
        index_root = find(rootTooth.vertices(:,2)<=(centerRoot(2) - (crownYdir(1)*(rootTooth.vertices(:,1) ...
                         - centerRoot(1))+crownYdir(3)*(rootTooth.vertices(:,3) ...
                         - centerRoot(3)))/crownYdir(2)));       % 牙齿自身坐标系下，y轴高度高于重心的点。
        vz_tooth = abs(centerRoot(2) - (crownYdir(1)*(rootTooth.vertices(index_root,1)-centerRoot(1))...
                     +crownYdir(3)*(rootTooth.vertices(index_root,3)-centerRoot(3)))/crownYdir(2));
        in_tooth = find(vz_tooth == max(vz_tooth));
        point_tooth = rootTooth.vertices(index_root(in_tooth),:);
        centcrownintooth = centerCrown + point_tooth - point_crown;  % 切后标准根的边缘中心
        
        % for debug
        OBJwriteVertices('centerRoot.obj', centerRoot);  
        OBJwriteVertices('point_tooth.obj', point_tooth);  
        OBJwriteVertices('切后标准根的边缘中心.obj', centcrownintooth);  
        
        
       % 2. 对齐——利用标准牙和病人牙齿的中心点以及三轴进行对齐，对齐后按照病人牙冠的切割位置进行切割
        R = inv([crown_ax.x; crown_ax.y; crown_ax.z]) * axi;   % 旋转
        C = centcrownintooth * R;         
        tooth_T = (rootTooth.vertices) *R + repmat((centerPatient - C),length(rootTooth.vertices),1);
        
        %       2.1 旋转平移带根标准牙
        index_1 = find(patientTooth.vertex(:,2)>=(centerPatient(2) - (toothYdir(1)*(patientTooth.vertex(:,1) - p_bingren(1))+toothYdir(3)*(patientTooth.vertex(:,3) - p_bingren(3)))/toothYdir(2))...
                 & patientTooth.vertex(:,2)<=(p_bingren(2) - (toothYdir(1)*(patientTooth.vertex(:,1) - ...
                 p_bingren(1))+toothYdir(3)*(patientTooth.vertex(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        patientMergeVers = patientTooth.vertex(index_1,:);
        
        index_2 = find(tooth_T(:,2)>=(centerPatient(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2))...
                & tooth_T(:,2)<=(p_bingren(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        rootMergeVers = tooth_T(index_2,:);
        
        % for debug
         OBJwriteVertices('2.1旋转平移后的带根标准牙.obj', tooth_T);
         OBJwriteVertices('2.1病人牙冠上融合区域的点.obj', patientMergeVers);
         OBJwriteVertices('2.1带根标准牙融合区域的点.obj', rootMergeVers);
        
        %       2.2 确定带根标准牙的切割部分。
        for i = 1:length(rootMergeVers)
           [minValue,r]=mindis(patientMergeVers, rootMergeVers(i,:),1);
           minvalue(i) = minValue;  row(i) = r;
        end
        point11 = patientMergeVers(row,:);
        [R,t,BRt,e,~,~] = icp(point11,rootMergeVers);
        tooth_T = bsxfun(@plus,tooth_T*R,t);        
        index_tooth =  find(tooth_T(:,2)>(p_bingren(2) - ( toothYdir(1)*(tooth_T(:,1) ...
                         - p_bingren(1))+ toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/ toothYdir(2))-2);  %切出来的牙根的索引值%注意！按照病人牙冠方向切
        rootCutVers = tooth_T(index_tooth ,:);
      

        %       2.3 牙冠变形控制部分
        index_crown_change = find(patientTooth.vertex(:,2)>(centerPatient(2) - (toothYdir(1)*(patientTooth.vertex(:,1) - p_bingren(1))+toothYdir(3)*(patientTooth.vertex(:,3) - p_bingren(3)))/toothYdir(2))-1 ...
                                & patientTooth.vertex(:,2)<=(p_bingren(2) - (toothYdir(1)*(patientTooth.vertex(:,1) - p_bingren(1))+toothYdir(3)*(patientTooth.vertex(:,3) - p_bingren(3)))/toothYdir(2))+2);
        crownforchange = patientTooth.vertex(index_crown_change,:);
        
       % for debug
        OBJwriteVertices('切出来的牙根顶点.obj', rootCutVers);

        %       2.5 牙冠和牙根合并成一个网格网格？？？？
        mergedToothVers = [cutPatientVers; rootCutVers];
        [mergedToothFace]=MyCrustOpen(mergedToothVers);
        % for debug
        writeOBJ('合并网格.Obj', mergedToothVers, mergedToothFace);             % ！！！ 三角片朝向混乱

       
        %       2.6 网格补洞
        mergedToothBdr = select_holes_and_boundary(mergedToothVers, mergedToothFace);
        ff = fill_mesh_holes(mergedToothVers, mergedToothFace, mergedToothBdr,'closed',200);
        ff = double(ff);
        index_tooth_change =  find(mergedToothVers(:,2)>(centerPatient(2) - ( toothYdir(1)*(mergedToothVers(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedToothVers(:,3) - p_bingren(3)))/ toothYdir(2))-1 ...
            &mergedToothVers(:,2)<=(p_bingren(2) - ( toothYdir(1)*(mergedToothVers(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedToothVers(:,3) - p_bingren(3)))/ toothYdir(2))+1);
        pchange = mergedToothVers(index_tooth_change,:);
        
        
        %       2.7 病人牙冠部分多截一部分出来，用以确定标准牙根变形后点的位置，找到tooth_T中离crownforchange中最近点的位置，并将这些点变为新的点，即牙冠点的位置
        mergedToothVers_copy = mergedToothVers;
        BZ1 = zeros(size(mergedToothVers_copy,1),3);
        indices = 1:size(mergedToothVers_copy,1);
        
        %       2.8 不变形的网格
        exterior = indices(mergedToothVers(:,2)<(p_bingren(2) - ( toothYdir(1)*(mergedToothVers(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedToothVers(:,3) - p_bingren(3)))/ toothYdir(2))-0.6...
            |mergedToothVers(:,2)>=(p_bingren(2) - ( toothYdir(1)*(mergedToothVers(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedToothVers(:,3) - p_bingren(3)))/ toothYdir(2))+3);
    end
else                             % 下颌牙
    % 。。。。。。。
end



[Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(mergedToothVers_copy,1), ff, exterior);


[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(mergedToothVers,ff, 'ext', 'voronoi', 'no_flatten',Omega,N0,N1);


%找到牙根部分跟牙冠变形部分最近的点，将牙根上的点拉至牙冠处
for i = 1:length(pchange)
   [minValue,r]=mindis(crownforchange,pchange(i,:),1);
   minvalue(i) = minValue;  row(i) = r;
   mergedToothVers_copy(index_tooth_change(i),:) = crownforchange(r,:);
end

bi_V = biharm_solve_with_factor( ...
    bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
    ff, mergedToothVers_copy, Omega, N0, N1, 'ext', 'no_flatten', BZ1, mergedToothVers);


figure()
trisurf(ff,bi_V(:,1),bi_V(:,2),bi_V(:,3),'facecolor','c','edgecolor','b')
% axis image
namestr3 = [num2str(x),'.','obj'];
writeOBJ(namestr3,bi_V,ff)

 
