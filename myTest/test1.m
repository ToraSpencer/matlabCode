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
        root = dentalwithtooth1(i);         % 有牙根的标准牙网格
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
writeOBJ('标准牙冠.obj', crown.model.vertex, crown.model.face);
writeOBJ('带牙根标准牙.obj', root.vertices, root.faces);        % ！！！带牙根的标准牙网格三角片是反的，是否有问题？
crown_center = mean(crown.model.vertex);
crown_ax_xline = getDirLine(crown_ax.x, crown_center);
crown_ax_yline = getDirLine(crown_ax.y, crown_center);
crown_ax_zline = getDirLine(crown_ax.z, crown_center);
OBJwriteVertices('标准牙x轴.obj', crown_ax_xline);
OBJwriteVertices('标准牙y轴.obj', crown_ax_yline);
OBJwriteVertices('标准牙z轴.obj', crown_ax_zline);
crown_ax_dirs = [crown_ax.x; crown_ax.y; crown_ax.z];
OBJwriteVertices('标准牙牙轴方向向量.obj', crown_ax_dirs);




if s ==1 || s ==2
    fdi = textread('FDIUpper__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
        disp( '该病例没有此牙齿 ');
        return;
       
    else
        namestr1 = ['toothUpper_',num2str(toothIdx-1),'.','obj'];
        bingren = Read_Obj(namestr1);                       % 病人的牙冠网格
        namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);                       % 牙龈线点云
        axis = ReadObj('AXISUpper_.obj');
        axi = axis( (3*(toothIdx-1) + 1) : 3*toothIdx, :);       % 当前牙齿的三个牙轴方向向量
        OBJwriteVertices('病人牙轴方向向量.obj', axi);
        OBJwriteVertices('病人牙龈线.obj', gumline);
        %       1.2 处理
        bingren_center = mean(bingren.vertex);     % 病人牙冠网格中心
        ymax = max(gumline(:,2));             % 牙龈线点云中y坐标最大值
        p_bingren =[bingren_center(1), ymax, bingren_center(3)];
        toothYdir = axi(2,:);                    % 当前病人牙齿的y轴方向
        bingren_vers = bingren.vertex;
        index_bingren = find(bingren_vers(:,2) <= (p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) ...
                - p_bingren(1)) + toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))-2);    % 超参数为2
        point_bingren = bingren_vers(index_bingren,:);      % 牙齿自身坐标系下，只保留y方向上高于p_bingren点2mm的点。
        
        disp(size(index_bingren));      % for debug
        
        
       %% FOR DEBUG
        figure(1);
        hold on;
        plot3(bingren.vertex(:,1), bingren.vertex(:,2), bingren.vertex(:,3), '.');
        plot3(gumline(:,1), gumline(:, 2), gumline(:, 3), '*');
        plot3(p_bingren(:,1), p_bingren(:, 2), p_bingren(:, 3), '*');
        hold off;
        
        OBJwriteVertices('p_bingren.obj', p_bingren);
       
        moved_p_bingren = p_bingren + 2 * toothYdir; 
        figure(2);
        hold on;
        plot3(point_bingren(:,1), point_bingren(:, 2), point_bingren(:, 3), '.');
        plot3(moved_p_bingren(:,1), moved_p_bingren(:, 2), moved_p_bingren(:, 3), '*');
        hold off;
        
        writeOBJ('病人牙冠.obj', bingren.vertex, bingren.face);
        OBJwriteVertices('moved_p_bingren.obj', moved_p_bingren);
        OBJwriteVertices('除去底部的病人牙冠点集.obj', point_bingren);
        
        bingren_xline = getDirLine(axi(1,:), bingren_center);
        bingren_yline = getDirLine(axi(2,:), bingren_center);
        bingren_zline = getDirLine(axi(3,:), bingren_center);
        OBJwriteVertices('病人牙x轴.obj', bingren_xline);
        OBJwriteVertices('病人牙y轴.obj', bingren_yline);
        OBJwriteVertices('病人牙z轴.obj', bingren_zline);
        
       %%
        
        %       标准牙冠
        face_crown_biaozhun = crown.model.face;
        vertex_crown_biaozhun = crown.model.vertex;
        crown_biaozhun_center = mean(vertex_crown_biaozhun);     % 标准牙冠重心
        biaozhunYdir = crown_ax.y;
        index_crown_biaozhun = find(vertex_crown_biaozhun(:,2)<=(crown_biaozhun_center(2) - (biaozhunYdir(1)*(vertex_crown_biaozhun(:,1) ...
                    - crown_biaozhun_center(1))+biaozhunYdir(3)*(vertex_crown_biaozhun(:,3) - crown_biaozhun_center(3)))/biaozhunYdir(2)));
        v_crown_biaozhun = abs(crown_biaozhun_center(2) - (biaozhunYdir(1)*(vertex_crown_biaozhun(index_crown_biaozhun,1)-crown_biaozhun_center(1))...
                     +biaozhunYdir(3)*(vertex_crown_biaozhun(index_crown_biaozhun,3)-crown_biaozhun_center(3)))/biaozhunYdir(2));
        in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
        point_crown =  vertex_crown_biaozhun(index_crown_biaozhun(in_crown),:);
        
       %%  for debug
        OBJwriteVertices('point_crown.obj', point_crown);
       %%
       
        % 标准牙根
        rootVers = root.vertices;
        root_center = mean(rootVers);
        index_root = find(rootVers(:,2)<=(root_center(2) - (biaozhunYdir(1)*(rootVers(:,1) ...
                         - root_center(1))+biaozhunYdir(3)*(rootVers(:,3) ...
                         - root_center(3)))/biaozhunYdir(2)));       % 牙齿自身坐标系下，y轴高度高于重心的点。
        vz_tooth = abs(root_center(2) - (biaozhunYdir(1)*(rootVers(index_root,1)-root_center(1))...
                     +biaozhunYdir(3)*(rootVers(index_root,3)-root_center(3)))/biaozhunYdir(2));
        in_tooth = find(vz_tooth == max(vz_tooth));
        point_tooth = rootVers(index_root(in_tooth),:);
        centcrownintooth = crown_biaozhun_center + point_tooth - point_crown;  % 切后标准根的边缘中心
        
        %% for debug
        chosenRootVers = rootVers(index_root, :);
        OBJwriteVertices('root_center.obj', root_center);  
        OBJwriteVertices('chosenRootVers.obj', chosenRootVers);      % 牙齿自身坐标系下，y轴高度高于重心的点。
        OBJwriteVertices('point_tooth.obj', point_tooth);  
        OBJwriteVertices('切后标准根的边缘中心.obj', centcrownintooth);  
        %%
        
       %% 2. 对齐——利用标准牙和病人牙齿的中心点以及三轴进行对齐，对齐后按照病人牙冠的切割位置进行切割
        R = inv([crown_ax.x; crown_ax.y; crown_ax.z]) * axi;   % 旋转
        C = centcrownintooth * R;         
        tooth_T = (rootVers) *R + repmat((bingren_center - C),length(rootVers),1);
        
        %       2.1 整体变化牙齿的形态
        index_1 = find(bingren_vers(:,2)>=(bingren_center(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))...
                 & bingren_vers(:,2)<=(p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) - ...
                 p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        point1 = bingren_vers(index_1,:);
        
        index_2 = find(tooth_T(:,2)>=(bingren_center(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2))...
                & tooth_T(:,2)<=(p_bingren(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        point2 = tooth_T(index_2,:);
        
        %% for debug
         OBJwriteVertices('旋转平移后的带根标准牙.obj', tooth_T);
         OBJwriteVertices('point1.obj', point1);
         OBJwriteVertices('point2.obj', point2);
        
       %%
        %       2.2 确定变形参数
        for i = 1:length(point2)
           [minValue,r]=mindis(point1,point2(i,:),1);
           minvalue(i) = minValue;  row(i) = r;
        end
        point11 = point1(row,:);
        [R,t,BRt,e,~,~] = icp(point11,point2);
        tooth_T = bsxfun(@plus,tooth_T*R,t);        
        index_tooth =  find(tooth_T(:,2)>(p_bingren(2) - ( toothYdir(1)*(tooth_T(:,1) ...
                         - p_bingren(1))+ toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/ toothYdir(2))-2);%切出来的牙根的索引值%注意！按照病人牙冠方向切
        rootCutVers = tooth_T(index_tooth ,:);
        

        

        %       2.3 牙冠变形控制部分
        index_crown_change = find(bingren_vers(:,2)>(bingren_center(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))-1 ...
                                & bingren_vers(:,2)<=(p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))+2);
        crownforchange = bingren_vers(index_crown_change,:);
        
       %% for debug
        OBJwriteVertices('切出来的牙根顶点.obj', rootCutVers);
        figure()
        trimesh(MyCrustOpen(point_bingren),point_bingren(:,1),point_bingren(:,2),point_bingren(:,3))
        hold on
        trimesh(MyCrustOpen(rootCutVers),rootCutVers(:,1),rootCutVers(:,2),rootCutVers(:,3))
        hold off
       %%
        
        %       2.5 牙冠和牙根合并成一个网格网格？？？？
        mergedTooth = [point_bingren; rootCutVers];
        [mergedTooth_face]=MyCrustOpen(mergedTooth);
        figure()
        trisurf(mergedTooth_face,mergedTooth(:,1),mergedTooth(:,2),mergedTooth(:,3),'facecolor','c','edgecolor','b');
        
        %% for debug
        writeOBJ('合并网格.Obj', mergedTooth, mergedTooth_face);             % ！！！ 三角片朝向混乱
        %%
       
        %       2.6 网格补洞
        b = select_holes_and_boundary(mergedTooth, mergedTooth_face);
        ff = fill_mesh_holes(mergedTooth, mergedTooth_face, b,'closed',200);
        ff = double(ff);
        index_tooth_change =  find(mergedTooth(:,2)>(bingren_center(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))-1 ...
            &mergedTooth(:,2)<=(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))+1);
        pchange = mergedTooth(index_tooth_change,:);
        
        
        %       2.7 病人牙冠部分多截一部分出来，用以确定标准牙根变形后点的位置，找到tooth_T中离crownforchange中最近点的位置，并将这些点变为新的点，即牙冠点的位置
        p_T = mergedTooth;
        bi_V = zeros(size(p_T));
        bi_bndtype = 'ext';
        BZ1 = zeros(size(p_T,1),3);
        reduction = 'no_flatten';masstype = 'voronoi';
        indices = 1:size(p_T,1);
        
        %       2.8 不变形的网格
        exterior = indices(mergedTooth(:,2)<(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))-0.6...
            |mergedTooth(:,2)>=(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))+3);
    end
else
    fdi = textread('FDILower__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
       disp( '该病例没有此牙齿 ');
       return;
    else
        namestr1 = ['toothLower_',num2str(toothIdx-1),'.','obj'];
        bingren = Read_Obj(namestr1);
        namestr2 = ['gumlineLower_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);
        axis = ReadObj('AXISLower_.obj');
        axi = axis(3*(toothIdx-1)+1:3*toothIdx,:);
    
         %处理
        bingren_center = mean(bingren.vertex);
        y = min(gumline(:,2));
        p_bingren =[bingren_center(1),y(1),bingren_center(3)];toothYdir = axi(2,:);bingren_vers = bingren.vertex;
        index_bingren = find(bingren_vers(:,2)>=(p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))+1);
        point_bingren = bingren_vers(index_bingren,:);
        
        %标准牙冠
        face_crown_biaozhun = crown.model.face;vertex_crown_biaozhun = crown.model.vertex;
        crown_biaozhun_center = mean(vertex_crown_biaozhun);
        biaozhunYdir = crown_ax.y;
        index_crown_biaozhun = find(vertex_crown_biaozhun(:,2)>=(crown_biaozhun_center(2) - (biaozhunYdir(1)*(vertex_crown_biaozhun(:,1) ...
                    - crown_biaozhun_center(1))+biaozhunYdir(3)*(vertex_crown_biaozhun(:,3) - crown_biaozhun_center(3)))/biaozhunYdir(2)));
        for i =1:length(index_crown_biaozhun)
            v_crown_biaozhun(i) = abs(dot((vertex_crown_biaozhun(index_crown_biaozhun(i),:)-crown_biaozhun_center),biaozhunYdir)/sqrt(sum(biaozhunYdir.*biaozhunYdir)));
        end
        in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
        
        point_crown =  vertex_crown_biaozhun(index_crown_biaozhun(in_crown),:);
 

        %标准牙根
        rootVers = root.vertices;
        f = root.faces;
        root_center = mean(rootVers);
        index_root = find(rootVers(:,2)>=(root_center(2) - (biaozhunYdir(1)*(rootVers(:,1) ...
                         - root_center(1))+biaozhunYdir(3)*(rootVers(:,3) - root_center(3)))/biaozhunYdir(2)));
        for i =1:length(index_root)
            vz_tooth(i) = abs(dot((rootVers(index_root(i),:)-root_center),biaozhunYdir)/sqrt(sum(biaozhunYdir.*biaozhunYdir)));
        end
        in_tooth = find(vz_tooth == max(vz_tooth));
        
        point_tooth = rootVers(index_root(in_tooth),:);
        rootVers = rootVers-point_tooth+point_crown;
        
        
        %%利用标准牙和病人牙齿的中心点以及三轴进行对齐，对齐后按照病人牙冠的切割位置进行切割
        R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*axi;
        C=crown_biaozhun_center*R;
        tooth_T = (rootVers) *R + repmat((bingren_center - C),length(rootVers),1);
        figure()
        trisurf(f,tooth_T(:,1),tooth_T(:,2),tooth_T(:,3),'facecolor','c','edgecolor','b')
        hold on
        trisurf(bingren.face,bingren_vers(:,1),bingren_vers(:,2),bingren_vers(:,3),'facecolor','c','edgecolor','r')
        
        %整体变化牙齿的形态
        index_1 = find(bingren_vers(:,2)<=(bingren_center(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2)) ...
                & bingren_vers(:,2)>=(p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        point1 = bingren_vers(index_1,:);
        index_2 = find(tooth_T(:,2)<=(bingren_center(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2)) ...
                & tooth_T(:,2)>=(p_bingren(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        point2 = tooth_T(index_2,:);
        
        
        %确定变形参数
        for i = 1:length(point2)
           [minValue,r]=mindis(point1,point2(i,:),1);
           minvalue(i) = minValue;  row(i) = r;
        end
        point11 = point1(row,:);
        [R,t,BRt,e,~,~] = icp(point11,point2);
        tooth_T = bsxfun(@plus,tooth_T,t);
        
        index_tooth =  find(tooth_T(:,2)<(p_bingren(2) - ( toothYdir(1)*(tooth_T(:,1) ...
                         - p_bingren(1))+ toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/ toothYdir(2))+1);%切出来的牙根的索引值%注意！按照病人牙冠方向切
        rootCutVers = tooth_T(index_tooth ,:);
        
 
        %牙冠变形控制部分
        index_crown_change = find(bingren_vers(:,2)<(bingren_center(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))+1 ...
                                & bingren_vers(:,2)>=(p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))-1);
        crownforchange = bingren_vers(index_crown_change,:);
        
        %牙冠牙根分别形成网格
        [tooth_t]=MyCrustOpen(point_bingren);
        figure()
        trimesh(tooth_t,point_bingren(:,1),point_bingren(:,2),point_bingren(:,3))
        hold on
        [tooth_t2]=MyCrustOpen(rootCutVers); 
        trimesh(tooth_t2,rootCutVers(:,1),rootCutVers(:,2),rootCutVers(:,3))
        hold off
       
        %牙冠和牙根形成网格
        mergedTooth = [point_bingren;rootCutVers];
        [t]=MyCrustOpen(mergedTooth);
        figure()
        trisurf(t,mergedTooth(:,1),mergedTooth(:,2),mergedTooth(:,3),'facecolor','c','edgecolor','b')
        
        %网格补洞
        b = select_holes_and_boundary(mergedTooth,t);   % ！！！当前运行到这里会报错，可能是matlab2015不兼容的问题。
        ff = fill_mesh_holes(mergedTooth,t,b,'closed',200);
        ff = double(ff);
        index_tooth_change =  find(mergedTooth(:,2)<(bingren_center(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))+1 ...
            &mergedTooth(:,2)>=(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))-1);
        pchange = mergedTooth(index_tooth_change,:);
        
        %病人牙冠部分多截一部分出来，用以确定标准牙根变形后点的位置，找到tooth_T中离crownforchange中最近点的位置，并将这些点变为新的点，即牙冠点的位置
        p_T = mergedTooth;
        bi_V = zeros(size(p_T));
        bi_bndtype = 'ext';
        BZ1 = zeros(size(p_T,1),3);
        reduction = 'no_flatten';masstype = 'voronoi';
        indices = 1:size(p_T,1);
        
        %不变形的网格
        exterior = indices(mergedTooth(:,2)>(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))+1 ...
            |mergedTooth(:,2)<=(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))-3);
    end
    
    
    
end


[Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(p_T,1), ff, exterior);
[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(mergedTooth,ff, bi_bndtype,masstype,reduction,Omega,N0,N1);
%找到牙根部分跟牙冠变形部分最近的点，将牙根上的点拉至牙冠处
for i = 1:length(pchange)
   [minValue,r]=mindis(crownforchange,pchange(i,:),1);
   minvalue(i) = minValue;  row(i) = r;
   p_T(index_tooth_change(i),:) = crownforchange(r,:);
end
bi_V = biharm_solve_with_factor( ...
    bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
    ff, p_T, Omega, N0, N1, bi_bndtype, reduction,BZ1,mergedTooth);
figure()
trisurf(ff,bi_V(:,1),bi_V(:,2),bi_V(:,3),'facecolor','c','edgecolor','b')
% axis image
namestr3 = [num2str(x),'.','obj'];
writeOBJ(namestr3,bi_V,ff)

 
