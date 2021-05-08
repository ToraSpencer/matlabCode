clc
clear all
functionname='testforfusion.m'; functiondir=which(functionname);
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
% up = Read_Obj('117608_upperMesh.obj');%病人的牙齿
% upaxis = textread ('117608_upperDir_centre.txt');%病人的牙轴
% low = Read_Obj('117608_lowerMesh.obj');
% lowaxis = textread ('117608_lowerDir_centre.txt');
% for i = 1:length(up)
%    up(i).center = upaxis(4*(i-1)+1,:);
%    up(i).orientation = upaxis(4*(i-1)+2:4*i,:);
%    face = [];
%    l = length(up(i).vertex);
%    face = [l+up(i).face(:,3)+1,l+up(i).face(:,2)+1,l+up(i).face(:,1)+1];
%    up(i).face = face;
%    trimesh(face,up(i).vertex(:,1),up(i).vertex(:,3),up(i).vertex(:,2))
%    hold on
%    
% end 
% for i = 1:length(low)
%    low(i).center = lowaxis(4*(i-1)+1,:);
%    low(i).orientation = lowaxis(4*(i-1)+2:4*i,:);
%    face = [];
%    l = length(low(i).vertex);
%    face = [l+low(i).face(:,3)+1,l+low(i).face(:,2)+1,l+low(i).face(:,1)+1];
%    low(i).face = face;
%    trimesh(face,low(i).vertex(:,1),low(i).vertex(:,3),low(i).vertex(:,2))
%    hold on
%    
% end 
% % 加载标准牙
load('dental_crown.mat');
load('dentalmodelwithroot0.1.mat');
load('axisofdental_crowm');
load('trfromtoothtocrown.mat');
load('up_bingren.mat');
load('low_bingren.mat');

%相同的FDI号拿出来
prompt  = '请输入FDI编号 \n';
x = input(prompt);
for i = 1:28
    if  dentalwithtooth1(i).ID == x
        root = dentalwithtooth1(i);
    end
end
 
for j = 1:28
    if (upax(j).fid == x)
        crown_ax = upax(j);
    end
end
for k = 1:28
    if (dental_crown(k).fid == x)
        crown = dental_crown(k);
    end
end
s = fix(x/10);  %取整
g = mod(x,10);%取余
if s ==1 || s ==2
    if s == 2
       bingren = up(8+g); 
    end
    if s == 1
       bingren = up(9-g); 
    end
    %处理
    p_bingren = bingren.center;n_bingren = bingren.orientation(2,:);v_bingren = bingren.vertex;
    index_bingren = find(v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2)));
    point_bingren = v_bingren(index_bingren,:);
    %标准牙冠
    face_crown_biaozhun = crown.model.face;vertex_crown_biaozhun = crown.model.vertex;
    p_crown_biaozhun = mean(vertex_crown_biaozhun);
    n_crown_biaozhun = crown_ax.y;
    index_crown_biaozhun = find(vertex_crown_biaozhun(:,2)<=(p_crown_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_crown_biaozhun(:,1) ...
                - p_crown_biaozhun(1))+n_crown_biaozhun(3)*(vertex_crown_biaozhun(:,3) - p_crown_biaozhun(3)))/n_crown_biaozhun(2)));
    v_crown_biaozhun = abs(p_crown_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_crown_biaozhun(index_crown_biaozhun,1)-p_crown_biaozhun(1))...
                 +n_crown_biaozhun(3)*(vertex_crown_biaozhun(index_crown_biaozhun,3)-p_crown_biaozhun(3)))/n_crown_biaozhun(2));
    in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
    point_crown =  vertex_crown_biaozhun(index_crown_biaozhun(in_crown),:);
   
    %标准牙根
    vertex_tooth_biaozhun = root.vertices;
    f = root.faces;
    p_tooth_biaozhun = mean(vertex_tooth_biaozhun);
%     index_tooth_biaozhun = find(vertex_tooth_biaozhun(:,2)<=(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) ...
%                      - p_tooth_biaozhun(1))+n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - p_tooth_biaozhun(3)))/n_crown_biaozhun(2)));
    vz_tooth = abs(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1)-p_tooth_biaozhun(1))...
               +n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3)-p_tooth_biaozhun(3)))/n_crown_biaozhun(2));
    in_tooth = find(vz_tooth == max(vz_tooth));
    point_tooth = vertex_tooth_biaozhun(in_tooth,:);
    centcrownintooth = p_crown_biaozhun+point_tooth-point_crown;%切后标准根的边缘中心
    index_tooth =  find(vertex_tooth_biaozhun(:,2)>(centcrownintooth(2) - ( n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) ...
                     - centcrownintooth(1))+ n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - centcrownintooth(3)))/ n_crown_biaozhun(2)));%切出来的牙根的索引值%注意！按照病人牙冠方向切
%     plot3(vertex_tooth_biaozhun(index_tooth,1),vertex_tooth_biaozhun(index_tooth,3),vertex_tooth_biaozhun(index_tooth,2),'g*')
    
    tooth = vertex_tooth_biaozhun(index_tooth,:);
    R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*bingren.orientation;
    C=centcrownintooth*R;
    tooth_T = (tooth) *R + repmat((p_bingren - C),length(tooth),1);

    %牙冠变形控制部分
    index_crown_change = find(v_bingren(:,2)>(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-0.5...
                            & v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+2);
    crownforchange = v_bingren(index_crown_change,:);
    %牙冠牙根分别形成网格
    [tooth_t]=MyCrustOpen(point_bingren);
    figure()
    trimesh(tooth_t,point_bingren(:,1),point_bingren(:,2),point_bingren(:,3))
    hold on
    [tooth_t2]=MyCrustOpen(tooth); 
    trimesh(tooth_t2,tooth(:,1),tooth(:,2),tooth(:,3))
    hold off
    %牙冠和牙根形成网格
    p = [point_bingren;tooth_T];
    [t]=MyCrustOpen(p);
    figure()
    trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
    %网格补洞
    b = select_holes_and_boundary(p,t);
    ff = fill_mesh_holes(p,t,b,'closed',200);
    ff = double(ff);
    index_tooth_change =  find(p(:,2)>(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
        - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-0.5...
        &p(:,2)<=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
        - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+1);
    pchange = p(index_tooth_change,:);
    %病人牙冠部分多截一部分出来，用以确定标准牙根变形后点的位置，找到tooth_T中离crownforchange中最近点的位置，并将这些点变为新的点，即牙冠点的位置
    p_T = p;
    bi_V = zeros(size(p_T));
    bi_bndtype = 'ext';
    BZ1 = zeros(size(p_T,1),3);
    reduction = 'no_flatten';masstype = 'voronoi';
    indices = 1:size(p_T,1);
    %不变形的网格
    exterior = indices(p(:,2)<(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
        - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-0.6...
        |p(:,2)>=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
        - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+3);

else
    if s == 3
       bingren = low(8-g); 
    end
    if s == 4
       bingren = low(7+g); 
    end
     %处理
    p_bingren = bingren.center;n_bingren = bingren.orientation(2,:);v_bingren = bingren.vertex;
    index_bingren = find(v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2)));
    point_bingren = v_bingren(index_bingren,:);
    %标准牙冠
    face_crown_biaozhun = crown.model.face;vertex_crown_biaozhun = crown.model.vertex;
    p_crown_biaozhun = mean(vertex_crown_biaozhun);
    n_crown_biaozhun = crown_ax.y;
    index_crown_biaozhun = find(vertex_crown_biaozhun(:,2)>=(p_crown_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_crown_biaozhun(:,1) ...
                - p_crown_biaozhun(1))+n_crown_biaozhun(3)*(vertex_crown_biaozhun(:,3) - p_crown_biaozhun(3)))/n_crown_biaozhun(2)));
    v_crown_biaozhun = abs(p_crown_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_crown_biaozhun(index_crown_biaozhun,1)-p_crown_biaozhun(1))...
                 +n_crown_biaozhun(3)*(vertex_crown_biaozhun(index_crown_biaozhun,3)-p_crown_biaozhun(3)))/n_crown_biaozhun(2));
    in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
    point_crown =  vertex_crown_biaozhun(index_crown_biaozhun(in_crown),:);
    
    %标准牙根
    vertex_tooth_biaozhun = root.vertices;
    f = root.faces;
    p_tooth_biaozhun = mean(vertex_tooth_biaozhun);
%     index_tooth_biaozhun = find(vertex_tooth_biaozhun(:,2)>=(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) ...
%                      - p_tooth_biaozhun(1))+n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - p_tooth_biaozhun(3)))/n_crown_biaozhun(2)));
    vz_tooth = abs(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1)-p_tooth_biaozhun(1))...
               +n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3)-p_tooth_biaozhun(3)))/n_crown_biaozhun(2));
    in_tooth = find(vz_tooth == max(vz_tooth));
    point_tooth = vertex_tooth_biaozhun(in_tooth,:);
    centcrownintooth = p_crown_biaozhun+point_tooth -point_crown;%切后标准根的边缘中心
    index_tooth =  find(vertex_tooth_biaozhun(:,2)<(centcrownintooth(2) - ( n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) ...
                     - centcrownintooth(1))+ n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - centcrownintooth(3)))/ n_crown_biaozhun(2)));%切出来的牙根的索引值%注意！按照病人牙冠方向切
    %     plot3(vertex_tooth_biaozhun(index_tooth,1),vertex_tooth_biaozhun(index_tooth,3),vertex_tooth_biaozhun(index_tooth,2),'g*')
    tooth = vertex_tooth_biaozhun(index_tooth,:);
    
%     R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*bingren.orientation;
%     C=centcrownintooth*R;
%     tooth_T = (tooth) *R + repmat((p_bingren - C),length(tooth),1);

    %牙冠变形控制部分
    index_crown_change = find(v_bingren(:,2)<(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+0.5...
                            & v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-3);
    crownforchange = v_bingren(index_crown_change,:);
    index_tooth_change1 =  find(vertex_tooth_biaozhun(:,2)<(centcrownintooth(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) - centcrownintooth(1))+n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - centcrownintooth(3)))/n_crown_biaozhun(2))+0.5...
                            & vertex_tooth_biaozhun(:,2)>=(centcrownintooth(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) - centcrownintooth(1))+n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - centcrownintooth(3)))/n_crown_biaozhun(2))-2);
    toothforchange = vertex_tooth_biaozhun(index_tooth_change1,:);
    index_rand = randperm(length(crownforchange),length(toothforchange));%这边需要家一个判断
    crownforchange1= crownforchange(index_rand,:);
    [R,t,BRt,e,KDTA,KDTB] = icp(crownforchange1,toothforchange);
    tooth_T = bsxfun(@plus,tooth*R,t);
    %牙冠牙根分别形成网格
    [tooth_t]=MyCrustOpen(point_bingren);
    figure()
    trimesh(tooth_t,point_bingren(:,1),point_bingren(:,2),point_bingren(:,3))
    hold on
    [tooth_t2]=MyCrustOpen(tooth); 
    trimesh(tooth_t2,tooth(:,1),tooth(:,2),tooth(:,3))
    hold off
    %牙冠和牙根形成网格
    p = [point_bingren;tooth_T];
    [t]=MyCrustOpen(p);
    figure()
    trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
    %网格补洞
    b = select_holes_and_boundary(p,t);
    ff = fill_mesh_holes(p,t,b,'closed',200);
    ff = double(ff);
    index_tooth_change =  find(p(:,2)<(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
        - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+0.5...
        &p(:,2)>=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
        - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-3);
    pchange = p(index_tooth_change,:);
    %病人牙冠部分多截一部分出来，用以确定标准牙根变形后点的位置，找到tooth_T中离crownforchange中最近点的位置，并将这些点变为新的点，即牙冠点的位置
    p_T = p;
    bi_V = zeros(size(p_T));
    bi_bndtype = 'ext';
    BZ1 = zeros(size(p_T,1),3);
    reduction = 'no_flatten';masstype = 'voronoi';
    indices = 1:size(p_T,1);
    %不变形的网格
    exterior = indices(p(:,2)>(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
        - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+0.6...
        |p(:,2)<=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
        - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-3);
    
    
    
end


[Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(p_T,1), ff, exterior);
[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(p,ff, bi_bndtype,masstype,reduction,Omega,N0,N1);
%找到牙根部分跟牙冠变形部分最近的点，将牙根上的点拉至牙冠处
for i = 1:length(crownforchange)
   [minValue,r]=mindis(pchange,crownforchange(i,:),1);
   minvalue(i) = minValue;  row(i) = r;
   p_T(index_tooth_change(r),:) = crownforchange(i,:);
end
bi_V = biharm_solve_with_factor( ...
    bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
    ff, p_T, Omega, N0, N1, bi_bndtype, reduction,BZ1,p);
figure()
trisurf(ff,bi_V(:,1),bi_V(:,2),bi_V(:,3),'facecolor','c','edgecolor','b')



% for i = 1:14
%     figure()
%     trimesh(dentalwithtooth(i).model.face,dentalwithtooth(i).model.vertex(:,1),dentalwithtooth(i).model.vertex(:,2),dentalwithtooth(i).model.vertex(:,3))
%     hold on
%     trimesh(up(i+1).face,up(i+1).vertex(:,1),up(i+1).vertex(:,2),up(i+1).vertex(:,3))
%     %切出牙根
%     
%     
%     %切出病人牙冠,沿着长轴方向，在与长轴垂直的平面上，切的平面上的点确定问题，暂时是牙冠的中心点，至于牙根
%     
%     quiver3(up(i+1).center(1),up(i+1).center(2),up(i+1).center(3),up(i+1).orientation(3,1),up(i+1).orientation(3,2),up(i+1).orientation(3,3),10); 
%     
%     %牙轴对齐（包含一个固定点和轴，点是平面的中心点）
%     
% end

