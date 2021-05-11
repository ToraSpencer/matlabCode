clc
clear all
functionname='testfornewcrown.m'; functiondir=which(functionname);
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
% addpath([functiondir '病例2带有牙龈线数据\CM210127095\l'])
addpath([functiondir 'mesh process'])
addpath([functiondir '病例2带有牙龈线数据\CM210127095\u'])
addpath([functiondir 'modelStadard'])
addpath([functiondir 'testdata'])
load('dentalmodelwithroot0.1.mat')
load('dental_crown.mat')
% local = textread('location__.txt');%每颗牙的位置
fdi = textread('FDI__.txt');
prompt  = '请输入FDI编号 \n';
x = input(prompt);
n = find (fdi == x);
for i = 1:28
    if  dentalwithtooth1(i).ID == x
        root = dentalwithtooth1(i);
    end
end
for i = 1:28
    if (dental_crown(i).fid == x)
        crown = dental_crown(i);
    end
end
namestr1 = ['TOOTH_',num2str(n-1),'.','obj'];
tooth_crown = Read_Obj(namestr1);
namestr2 = ['gumline_',num2str(n-1),'.','obj'];
gumline =  ReadObj(namestr2);
axis = ReadObj('AXIS_.obj');
axi = axis(3*(n-1)+1:3*n,:);
%标准牙和病人的牙齿对齐
R = inv([root.nx;root.ny;root.nz])*axi;
C = mean(tooth_crown.vertex);
p_bingren = [C(1),max(gumline(:,2)),C(3)]; %gumline(find(gumline(:,2) == max(gumline(:,2))),:);
n_bingren = axi(2,:);
crown_biaozhun = (crown.model.vertex) *R;c=mean(crown_biaozhun);
crown_biaozhun = crown_biaozhun+C-c;
% % if length(crown_biaozhun)<=length(tooth_crown.vertex)
% %     [error,Reallignedsource,transform]=rigidICP(tooth_crown.vertex,crown_biaozhun,1,tooth_crown.face,crown.model.face);
% %     toothwithroot = transform.b*(root.vertices*R+C-c)*transform.T + transform.c(1,:);
% % else
% %     toothwithroot = root.vertices*R+C-c;
% % end


toothwithroot = root.vertices*R+C-c;


%标准牙变形
%控制变形位置为标准牙牙冠上方的位置，所以，先要找到牙冠（病人和标准牙的），将标准牙冠的点变形到病人的牙冠上，牙根位置保持不变
%% 因为标准牙模型的点比病人牙齿疏，所以，只要找离牙龈线最近的点，若是标准牙的点较密，则取最近的前几个点
%找牙冠上的点，用于判断牙冠跟牙根点
%找到标准牙冠上的点，并且确定牙根上不变的点
xmax_biaozhun = max(toothwithroot(:,1));xmin_biaozhun = min(toothwithroot(:,1));
xmid_biaozhun = mean([xmax_biaozhun;xmin_biaozhun]);
index_biaozhun_x = find(toothwithroot(:,1)>= xmid_biaozhun-0.2 & toothwithroot(:,1)<= xmid_biaozhun+0.2);
index_biaozhun = find(toothwithroot(index_biaozhun_x,2) == min(toothwithroot(index_biaozhun_x,2)));
row = [];k =3;
for i = 1:length(gumline)
   [~,r]=mindis(toothwithroot,gumline(i,:),k);%k表示排行k个最小
   row = [row;r];
end
%牙冠上的点与牙根上的点分开
[in_biaozhun] = segmentation_region_grow(toothwithroot,root.faces,toothwithroot(row,:));
if isempty(find(root.faces(in_biaozhun,:)==index_biaozhun_x(index_biaozhun)))
    kk_tooth_biaozhun = in_biaozhun;%kk_biaozhun牙根上的面片
    kk1 = (1:length(root.faces))';
    kk_crown_biaozhun  = setdiff(kk1, in_biaozhun);
else
    kk1 = (1:length(root.faces))';
    kk_tooth_biaozhun  = setdiff(kk1, in_biaozhun);
    kk_crown_biaozhun = in_biaozhun;
end



%提取病人牙冠
xmax_crown = max(tooth_crown.vertex(:,1));xmin_crown = min(tooth_crown.vertex(:,1));
xmid_crownn = mean([xmax_crown;xmin_crown]);
index_crown_x = find(tooth_crown.vertex(:,1)>= xmid_crownn-0.2 & tooth_crown.vertex(:,1)<= xmid_crownn+0.2);
index_crown = find(tooth_crown.vertex(index_crown_x,2) == min(tooth_crown.vertex(index_crown_x,2)));%取最大值和最小值取决于牙冠的位置

[in_bingren] = segmentation_region_grow(tooth_crown.vertex,tooth_crown.face,gumline);
if isempty(find(tooth_crown.face(in_bingren,:)==index_crown_x(index_crown)))
    kk1 = (1:length(tooth_crown.face))';
    kk_bingren  = setdiff(kk1, in_bingren);%这边需要重新确认
else
    kk_bingren = in_bingren;
    
end
%这边变形是否可以考虑牙根的边界点与牙冠边界点的对齐
%判断三角形一边的位置
tri = root.faces(kk_crown_biaozhun,:);
boundaries = detect_mesh_holes_and_boundary(tri);
innu = unique(reshape(boundaries{1},[],1));
point1 = toothwithroot(innu,:);
% lin = [tri(:,1:2);tri(:,2:3);[tri(:,1) tri(:,3)]];
% lin1 = sortrows(lin);
% [dataUnique,r]=unique(lin1,'rows');
% [dataUnique,r1]=unique(lin1,'rows','last');
% nu = r1 -r+1;
% innu = unique(reshape(dataUnique(nu == 1,:),[],1));
% [li,IA,IC] = unique(lin,'rows');
% B=tabulate(IC);
% nu = B((B(:,2)==1),1);
% innu = unique(reshape(li(nu,:),[],1));
point1 = toothwithroot(innu,:);
figure()
trisurf(root.faces,toothwithroot(:,1),toothwithroot(:,2),toothwithroot(:,3),'facecolor','c','edgecolor','b')
hold on
plot3(point1(:,1),point1(:,2),point1(:,3),'r*')
axis image


% in_crown_biaozhun = unique(reshape(root.faces(kk_crown_biaozhun,:),[],1));
indices = 1:size(toothwithroot,1);
% % in_crown_biaozhun = indices(toothwithroot(:,2)<(p_bingren(2) - ( n_bingren(1)*(toothwithroot(:,1) ...
% %             - p_bingren(1))+ n_bingren(3)*(toothwithroot(:,3) - p_bingren(3)))/ n_bingren(2))& toothwithroot(:,2)>(p_bingren(2) - ( n_bingren(1)*(toothwithroot(:,1) ...
% %             - p_bingren(1))+ n_bingren(3)*(toothwithroot(:,3) - p_bingren(3)))/ n_bingren(2)));
% point1 = toothwithroot(in_crown_biaozhun,:);
point1 = toothwithroot(row,:);
% % in = 1:size(tooth_crown.vertex,1);
% % in_crown_bingren  = in(tooth_crown.vertex(:,2)<(p_bingren(2) - ( n_bingren(1)*(tooth_crown.vertex(:,1) ...
% %             - p_bingren(1))+ n_bingren(3)*(tooth_crown.vertex(:,3) - p_bingren(3)))/ n_bingren(2))+0.5& tooth_crown.vertex(:,2)>(p_bingren(2) - ( n_bingren(1)*(tooth_crown.vertex(:,1) ...
% %             - p_bingren(1))+ n_bingren(3)*(tooth_crown.vertex(:,3) - p_bingren(3)))/ n_bingren(2)));
% % point2 = tooth_crown.vertex(in_crown_bingren,:);
in_crown_bingren = unique(reshape(tooth_crown.face(kk_bingren,:),[],1));
point2 = tooth_crown.vertex(in_crown_bingren,:);
tooth = toothwithroot;
for i = 1:length(point1)
    %找到与病人长轴方向垂直的平面
%     B =point2 -point1(i,:);
%     for j = 1:length(point2)
%         K(j) = dot(n_bingren,B(j,:))/sqrt(sum(dot(n_bingren,n_bingren))*sum(dot(B(j,:),B(j,:))));
%     end
%     r = find(K >= 0);
%     n = find(K(r) == min(K(r)));
    [~,r]=mindis(point2,point1(i,:),1);
    tooth(row(i),:) = point2(r,:);
% %     tooth(in_crown_biaozhun(i),:) = point2(r,:);
end
bi_V = zeros(size(toothwithroot));
bi_bndtype = 'ext';
BZ1 = zeros(size(toothwithroot,1),3);
reduction = 'no_flatten';masstype = 'voronoi';

% exterior = unique(reshape(root.faces(kk_unchange,:),[],1));

exterior = indices(toothwithroot(:,2)>(p_bingren(2) - ( n_bingren(1)*(toothwithroot(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(toothwithroot(:,3) - p_bingren(3)))/ n_bingren(2))+2 | ...
            toothwithroot(:,2)<(p_bingren(2) - ( n_bingren(1)*(toothwithroot(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(toothwithroot(:,3) - p_bingren(3)))/ n_bingren(2))-1);
[Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(toothwithroot,1), root.faces, exterior);
[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(toothwithroot,root.faces, bi_bndtype,masstype,reduction,Omega,N0,N1);
bi_V = biharm_solve_with_factor( bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M,root.faces, tooth, Omega, N0, N1, bi_bndtype, reduction,BZ1,toothwithroot);
figure()
trisurf(root.faces,bi_V(:,1),bi_V(:,2),bi_V(:,3),'facecolor','c','edgecolor','b')
hold on
trisurf(tooth_crown.face,tooth_crown.vertex(:,1),tooth_crown.vertex(:,2),tooth_crown.vertex(:,3),'facecolor','y','edgecolor','g')
plot3(gumline(:,1),gumline(:,2),gumline(:,3),'r*')
plot3(toothwithroot(row,1),toothwithroot(row,2),toothwithroot(row,3),'y*')
hold off
axis image
row = [];k =4;
for i = 1:length(gumline)
   [~,r]=mindis(bi_V,gumline(i,:),k);%k表示排行k个最小
   row = [row;r];
end
%牙冠上的点与牙根上的点分开
[in_biao] = segmentation_region_grow(bi_V,root.faces,bi_V(row,:));
if isempty(find(root.faces(in_biao,:)==index_biaozhun_x(index_biaozhun)))
    kk = in_biao;%kk_biaozhun牙根上的面片   
else
    kk1 = (1:length(root.faces))';
    kk  = setdiff(kk1, in_biao);   
end

tooth_e = bi_V(unique(reshape(root.faces(kk,:),[],1)),:);
crown_e = tooth_crown.vertex(unique(reshape(tooth_crown.face(kk_bingren,:),[],1)),:);

figure()
plot3(tooth_e(:,1),tooth_e(:,2),tooth_e(:,3),'r*')
hold on
plot3(crown_e(:,1),crown_e(:,2),crown_e(:,3),'g*')
plot3(bi_V(row,1),bi_V(row,2),bi_V(row,3),'y*')
hold off
axis image

p = [crown_e;tooth_e];
[t]=MyCrustOpen(p);
figure()
trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
b = select_holes_and_boundary(p,t);
ff = fill_mesh_holes(p,t,b,'closed',200);
figure()
trisurf(ff,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
axis image












% n_crown = axi(2,:);n_biaozhun = root.ny;
% p_crown = mean(tooth_crown.vertex);p_biaozhun = mean(root.vertices);
% %
% index_crown = find(tooth_crown.vertex(:,2)<=(p_crown(2)));
% v_crown = abs(p_crown(2) - (n_crown(1)*(tooth_crown.vertex(index_crown,1)-p_crown(1))+n_crown(3)...
%            *(tooth_crown.vertex(index_crown,3)-p_crown(3)))/n_crown(2));
% in_crown = find(v_crown == max(v_crown));
% point_crown =  tooth_crown.vertex(index_crown(in_crown),:);
% %
% index_biaozhun = find(root.vertices(:,2)<=(p_biaozhun(2)));
% v_biaozhun = abs(p_biaozhun(2) - (n_biaozhun(1)*(root.vertices(index_biaozhun,1)-p_biaozhun(1))+n_biaozhun(3)...
%            *(root.vertices(index_biaozhun,3)-p_biaozhun(3)))/n_biaozhun(2));
% in_biaozhun = find(v_biaozhun == max(v_biaozhun));
% point_biaozhun =  root.vertices(index_biaozhun(in_biaozhun),:);
% %平移
% root_v = root.vertices+point_crown-point_biaozhun;
%% 判断那个长那个短
if length(tooth_crown.vertex) > length(crown.model.vertex)
    index_rand = randperm(length(tooth_crown.vertex),length(crown.model.vertex));%这边需要家一个判断
    v = tooth_crown.vertex(index_rand,:);
end
    
[R,t] = icp(v,crown.model.vertex);
tooth_C = bsxfun(@plus,root.vertices*R,t);

%% 唇舌侧和近远中方向变化
%问题：以哪个点为中心确定厚度（每个不同位置，厚度不同）
%沿着唇舌侧和近远中方向，怎么定义这两个的边界
%考虑将牙龈线投影到二维上，在二维上发生变化



R = inv([root.nx;root.ny;root.nz])*axi;
tooth_T = ( root.vertices) *R;
xmax_crown = max(tooth_crown.vertex(:,1));xmin_crown = min(tooth_crown.vertex(:,1));
xmid_crownn = mean([xmax_crown;xmin_crown]);
index_crown_x = find(tooth_crown.vertex(:,1)>= xmid_crownn-0.2 & tooth_crown.vertex(:,1)<= xmid_crownn+0.2);
index_crown = find(tooth_crown.vertex(index_crown_x,2) == min(tooth_crown.vertex(index_crown_x,2)));%取最大值和最小值取决于牙冠的位置
point_crown = tooth_crown.vertex(index_crown_x(index_crown),:);

xmax_biaozhun = max(tooth_T(:,1));xmin_biaozhun = min(tooth_T(:,1));
xmid_biaozhun = mean([xmax_biaozhun;xmin_biaozhun]);
index_biaozhun_x = find(tooth_T(:,1)>= xmid_biaozhun-0.2 & tooth_T(:,1)<= xmid_biaozhun+0.2);
index_biaozhun = find(tooth_T(index_biaozhun_x,2) == min(tooth_T(index_biaozhun_x,2)));
point_biaozhun = tooth_T(index_biaozhun_x(index_biaozhun),:);
root_v = tooth_T+point_crown-point_biaozhun;
% figure()
% trimesh(tooth_crown.face,tooth_crown.vertex(:,1),tooth_crown.vertex(:,2),tooth_crown.vertex(:,3))
% hold on
% trimesh(root.faces,root_v(:,1),root_v(:,2),root_v(:,3))




for i = 1:length(fdi)
    seed = [];loc =[];
    namestr1 = ['TOOTH_',num2str(i-1),'.','obj'];
    tooth(i) = Read_Obj(namestr1);
    namestr2 = ['gumline_',num2str(i-1),'.','obj'];
    gumline{i} =  ReadObj(namestr2);
    [in] = segmentation_region_grow(tooth(i).vertex,tooth(i).face,gumline{i});
    
    seed = tooth(i).face(in,:);
    loc = unique([seed(:,1);seed(:,2);seed(:,3)]);
    crownpointslocal{i} = loc;
    axi{i} = axis(3*(i-1)+1:3*i,:);
    
    
    
    
    
    
end
    