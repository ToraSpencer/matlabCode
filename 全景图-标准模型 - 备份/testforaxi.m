% clc
% clear all
functionname='testforaxi.m'; functiondir=which(functionname);
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
up = Read_Obj('117608_upperMesh.obj');%病人的牙齿
addpath([functiondir 'modelStadard'])
addpath([functiondir 'testdata'])
upaxis = textread ('117608_upperDir_centre.txt');%病人的牙轴
% for i = 1:length(up)
%    up(i).center = upaxis(4*(i-1)+1,:);
%    up(i).orientation = upaxis(4*(i-1)+2:4*i,:);
%    face = [];
%    l = length(up(i).vertex);
%    face = [l+up(i).face(:,3)+1,l+up(i).face(:,2)+1,l+up(i).face(:,1)+1];
%    up(i).face = face;
%    trimesh(face,up(i).vertex(:,1),up(i).vertex(:,3),up(i).vertex(:,2))
%    hold on  
% end  
% % % 加载标准牙
% load('dental_crown.mat');
% load('dentalwithtooth.mat');
% load('axisofdental_crowm');
% load('trfromtoothtocrown.mat');
k = 1;
figure()
for g = 1:4
    for h = 1:7
         dental_crown(k).fid = 10*g+h;dentalwithtooth(k).fid = 10*g+h;
         namestr1 = ['tooth',num2str(g),num2str(h),'.','obj'];
         dental_crown(k).model = Read_Obj(namestr1);
         namestr2 = [num2str(g),'-',num2str(h),'.','stl'];
         dentalwithtooth0(k).model = stlread(namestr2);
         trimesh(dental_crown(k).model.face,dental_crown(k).model.vertex(:,1),dental_crown(k).model.vertex(:,2),dental_crown(k).model.vertex(:,3))
         hold on 
         trimesh(dentalwithtooth(k).model.face,dentalwithtooth(k).model.vertex(:,1),dentalwithtooth(k).model.vertex(:,2),dentalwithtooth(k).model.vertex(:,3))
         k = k+1;      
    end
end
% 分别找牙冠距离最小的点，求两点的位置，用来变换
% 牙根
pointpair = [];
for n = length(dentalwithtooth)/2+1:length(dentalwithtooth)-1
    %%通过遍历找最近的点
    for i = 1:length(dentalwithtooth(n).model.vertex)
        [minValue,r]=matchest(dentalwithtooth(n+1).model.vertex,dentalwithtooth(n).model.vertex(i,:));
        value(i) = minValue;row(i) = r;
    end
    k = find(value == min(value));m = row(k);
    pointpair = [pointpair;dentalwithtooth(n).model.vertex(k,:),dentalwithtooth(n+1).model.vertex(m,:)];
    pointfindedtooth(n-length(dentalwithtooth)/2,:) = mean([dentalwithtooth(n).model.vertex(k,:);dentalwithtooth(n+1).model.vertex(m,:)]); 
    plot3(pointfindedtooth(n-length(dentalwithtooth)/2,1),pointfindedtooth(n-length(dentalwithtooth)/2,2),pointfindedtooth(n-length(dentalwithtooth)/2,3), 'r*','MarkerSize',100) 
    hold on
    value = [];
    row = [];
end
hold on
% %牙冠
pointpair = [];
for n = length(dental_crown)/2+1:length(dental_crown)-1
    %%通过遍历找最近的点
    for i = 1:length(dental_crown(n).model.vertex)
        [minValue,r]=matchest(dental_crown(n+1).model.vertex,dental_crown(n).model.vertex(i,:));
        value(i) = minValue;row(i) = r;
    end
    k = find(value == min(value));m = row(k);
    pointpair = [pointpair;dental_crown(n).model.vertex(k,:),dental_crown(n+1).model.vertex(m,:)];
    pointfindedcr(n-length(dental_crown)/2,:) = mean([dental_crown(n).model.vertex(k,:);dental_crown(n+1).model.vertex(m,:)]); 
    plot3(pointfindedcr(n-length(dental_crown)/2,1),pointfindedcr(n-length(dental_crown)/2,2),pointfindedcr(n-length(dental_crown)/2,3), 'B*','MarkerSize',100) 
    value = [];
    row = [];
end
[d, z, tform] = procrustes(pointfindedcr, pointfindedtooth, 'Reflection',false);
% % houzhui = '.obj';
% % upaxis2 = textread ('toothaxeslog.txt','%s','delimiter','\n');
% % k =1;
% % for j = length(upaxis2)-167:6:length(upaxis2)
% %     strhead = deblank(upaxis2{j,1});Shead = regexp(strhead, '\s+', 'split');
% %     strx = deblank(upaxis2{j+3,1});stry = deblank(upaxis2{j+4,1});strz = deblank(upaxis2{j+5,1});
% %     Sx = regexp(strx, '\s+', 'split');Sy = regexp(stry, '\s+', 'split');Sz = regexp(strz, '\s+', 'split');
% %     upax(k).fid = str2num(Shead{1,3});upax(k).x = [str2num(Sx{1,3}),str2num(Sx{1,4}),str2num(Sx{1,5})];
% %     upax(k).y = [str2num(Sy{1,3}),str2num(Sy{1,4}),str2num(Sy{1,5})];upax(k).z = [str2num(Sz{1,3}),str2num(Sz{1,4}),str2num(Sz{1,5})];
% %     k = k+1;
% % end  
% % % 标准牙的变换
% % 
% for i = 1:14
%     dentalwithtooth(i).model.vertex = dentalwithtooth(i).model.vertex*tform.T+tform.c(1,:);
% end
for i = 15:28
    dentalwithtooth1(i).model.vertex = dentalwithtooth0(i).model.vertex*tform.T+tform.c(1,:);
    trimesh(dental_crown(i).model.face,dental_crown(i).model.vertex(:,1),dental_crown(i).model.vertex(:,2),dental_crown(i).model.vertex(:,3))
    hold on 
    trimesh(dentalwithtooth0(i).model.face,dentalwithtooth(i).model.vertex(:,1),dentalwithtooth(i).model.vertex(:,2),dentalwithtooth(i).model.vertex(:,3))
end
%%标准牙牙根与病人牙冠的融合，单颗牙做案例
for i = 5:length(up)
%     [~,MC]=curvatures(up(i).vertex(:,1),up(i).vertex(:,2),up(i).vertex(:,3),up(i).face);
    %病人牙冠
    index_bingren = [];v_bingren = [];
    v_bingren = up(i).vertex;
%     figure()
%     trimesh(up(i).face,up(i).vertex(:,1),up(i).vertex(:,3),up(i).vertex(:,2))
%     hold on
%     quiver3(p(1),p(3),p(2),n(1),n(3),n(2),10)
    p_bingren = up(i).center;
    n_bingren = up(i).orientation(2,:);
    y = 3:0.5:15;
    x = -24:0.5:-12;
%     [X,Y]=meshgrid(x,y);
%     Z = p(2) - (n(1)*(X-p(1))+n(3)*(Y-p(3)))/n(2);
%     mesh(X,Y,Z)
%     surf(X,Y,Z)
%     hold off
    index_bingren = find(v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2)));
    point_bingren = v_bingren(index_bingren,:);

%     figure()
%     plot3(point_bingren(:,1),point_bingren(:,3),point_bingren(:,2),'r*')
    %标准牙
    %牙冠部分
    face_crown_biaozhun = dental_crown(4).model.face;
    vertex_crown_biaozhun = dental_crown(4).model.vertex;
    p_crown_biaozhun = mean(vertex_crown_biaozhun);
    n_crown_biaozhun = upax(4).y;
    index_crown_biaozhun = find(vertex_crown_biaozhun(:,2)<=(p_crown_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_crown_biaozhun(:,1) ...
    - p_crown_biaozhun(1))+n_crown_biaozhun(3)*(vertex_crown_biaozhun(:,3) - p_crown_biaozhun(3)))/n_crown_biaozhun(2)));
    v_crown_biaozhun = abs(p_crown_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_crown_biaozhun(index_crown_biaozhun,1)-p_crown_biaozhun(1))...
                     +n_crown_biaozhun(3)*(vertex_crown_biaozhun(index_crown_biaozhun,3)-p_crown_biaozhun(3)))/n_crown_biaozhun(2));
    in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
    point_crown =  vertex_crown_biaozhun(index_crown_biaozhun(in_crown),:);
    
    %牙根部分
    vertex_tooth_biaozhun = dentalwithtooth(4).model.vertex;
    f = dentalwithtooth(4).model.face;
    p_tooth_biaozhun = mean(vertex_tooth_biaozhun);
    index_tooth_biaozhun = find(vertex_tooth_biaozhun(:,2)<=(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) ...
                     - p_tooth_biaozhun(1))+n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - p_tooth_biaozhun(3)))/n_crown_biaozhun(2)));
    vz_tooth = abs(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(index_tooth_biaozhun,1)-p_tooth_biaozhun(1))...
               +n_crown_biaozhun(3)*(vertex_tooth_biaozhun(index_tooth_biaozhun,3)-p_tooth_biaozhun(3)))/n_crown_biaozhun(2));
    in_tooth = find(vz_tooth == max(vz_tooth));
    point_tooth = vertex_tooth_biaozhun(index_tooth_biaozhun(in_tooth),:);
    centcrownintooth = p_crown_biaozhun+mean(point_tooth)-point_crown;%切后标准根的边缘中心
    index_tooth =  find(vertex_tooth_biaozhun(:,2)>(centcrownintooth(2) - ( n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) ...
                     - centcrownintooth(1))+ n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - centcrownintooth(3)))/ n_crown_biaozhun(2)));%切出来的牙根的索引值%注意！按照病人牙冠方向切
%     plot3(vertex_tooth_biaozhun(index_tooth,1),vertex_tooth_biaozhun(index_tooth,3),vertex_tooth_biaozhun(index_tooth,2),'g*')
    tooth = vertex_tooth_biaozhun(index_tooth,:);
    R = inv([upax(4).x;upax(4).y;upax(4).z])*up(5).orientation;
    C=centcrownintooth*R;
    tooth_T = (tooth) *R + repmat((p_bingren - C),length(tooth),1);
% %     tooth_crown =(point_bingren*[upax(4).x;upax(4).y;upax(4).z]-centcrownintooth+p_bingren)*inv(up(5).orientation);%这个是坐标系转换，而牙根和牙冠跟坐标系转换没有关系   
%     
% %     %移动牙根，进行数据对齐
% %     R = inv([upax(4).x;upax(4).y;upax(4).z])*up(5).orientation;
% % %     crown_bingrenT = point_bingren*R+centcrownintooth-p_bingren;
% %     C=centcrownintooth*R;
% % %     tooth_T = (tooth) *R - C+p_bingren;%此为2019版本
% %     tooth_T = (tooth) *R + repmat((p_bingren - C),length(tooth),1);
% %     plot3(tooth_T(:,1),tooth_T(:,3),tooth_T(:,2),'g*')
% %     %牙根部分和牙冠部分边缘重合
    index_crown_change = find(v_bingren(:,2)>(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-0.5...
                            & v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+1);
    crownforchange = v_bingren(index_crown_change,:);
%     index_tooth_change =  find(tooth_T(:,2)>(p_bingren(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2))...
%                             & tooth_T(:,2)<=(p_bingren(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2))+0.6);
%     toothforchange = vertex_tooth_biaozhun(index_tooth_change,:);
    
%     [R,t,BRt,e,KDTA,KDTB] = icp(crownforchange,toothforchange);
%     tooth_T_T = bsxfun(@plus,tooth*R,t);
%     index_rand = randperm(length(toothforchange),length(crownforchange));%这边需要家一个判断
%     randtoothforchang= toothforchange(index_rand,:);
% %     [~, ~, T] = procrustes(crownforchange,randtoothforchange, 'Reflection',false);
%     
    [tooth_t]=MyCrustOpen(point_bingren);
    figure()
    trimesh(tooth_t,point_bingren(:,1),point_bingren(:,2),point_bingren(:,3))
    hold on
    [tooth_t2]=MyCrustOpen(tooth); 
    trimesh(tooth_t2,tooth(:,1),tooth(:,2),tooth(:,3))
    hold off
    
    p = [point_bingren;tooth_T];
    [t]=MyCrustOpen(p);
    figure()
    trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
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
    exterior = indices(p(:,2)<(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
        - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-0.6...
        |p(:,2)>=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
        - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+2);
    [Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(p_T,1), ff, exterior);
    [bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(p,ff, bi_bndtype,masstype,reduction,Omega,N0,N1);
    for i = 1:length(crownforchange)
       [minValue,r]=mindis(pchange,crownforchange(i,:));
       minvalue(i) = minValue;  row(i) = r;
       p_T(index_tooth_change(r),:) = crownforchange(i,:);
    end
    bi_V = biharm_solve_with_factor( ...
        bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
        ff, p_T, Omega, N0, N1, bi_bndtype, reduction,BZ1,p);
    figure()
    trisurf(ff,bi_V(:,1),bi_V(:,2),bi_V(:,3),'facecolor','c','edgecolor','b')

    
    
    
    
    
    
    
    [iaface1,ibface1]=find(bsxfun(@eq,ff(:,1),index_tooth_change'));
    [iaface2,ibface2]=find(bsxfun(@eq,ff(iaface1,2),index_tooth_change'));
    [iaface3,ibface3]=find(bsxfun(@eq,ff(iaface1(iaface2),3),index_tooth_change'));
    yy = 1:length(ff);[iaF,ibF]=ismember(yy,iaface1(iaface2(iaface3)));
    iafacec = yy(iaF==0);
    [face_new,vertex_new,father] = remove_mesh_face(ff,p,iafacec);
    xyzn=meannorm_trismooth(vertex_new,face_new);
       

%% 对齐以后再进行边界检测
%     
%     [Indices_edgesS]=detectedges(sourceV,sourceF);
    %%边处理
%     F = up(i).face;index=index_bingren;
%     edges = [sort([ F(:,1)  F(:,2) ]')';sort([ F(:,1)  F(:,3) ]')';sort([ F(:,2)  F(:,3) ]')'];
%     [etemp1,ia,ic]=unique(edges,'rows','stable');
%     %去除切除了边
%     [ia1,ib1]=find(bsxfun(@eq,etemp1(:,1),index')); %边的第一个点的位置 ia:第一个向量中相同元素的位置 ib：对应的第二个想两种元素的位置
%     a1 = etemp1(ia1,2);
%     [ia2,ib2]=find(bsxfun(@eq,a1,index'));
%     %起点在去除点范围内，终点不在去除范围内，则考虑为边界点
%     xx = [1:length(ia1)]'; 
%     [ia,ib]=ismember(xx,ia2);
% %     [iax,ibx]=find(bsxfun(@eq,xx,ia2'));
%     pointscrown1 = etemp1(ia1(ia==0),2);
%     [ia3,ib3]=find(bsxfun(@eq,etemp1(:,2),index'));
%     a2 = etemp1(ia3,1);
%     [ia4,ib4]=find(bsxfun(@eq,a2,index'));
%     yy = [1:length(ia3)]';
%     [ia,ib]=ismember(yy,ia4);
%     pointscrown2 = etemp1(ia1(ia==0),1);
%     crownedgepoints = v_bingren(unique(pointscrown1),:);
%     %牙冠
%     [iaface1,ibface1]=find(bsxfun(@eq,F(:,1),index'));
%     [iaface2,ibface2]=find(bsxfun(@eq,F(iaface1,2),index'));
%     [iaface3,ibface3]=find(bsxfun(@eq,F(iaface1(iaface2),3),index'));
%    
%     yy = 1:length(F);[iaF,ibF]=ismember(yy,iaface1(iaface2(iaface3)));
%     iafacec = yy(iaF==0);
%     [face_new,vertex_new,father] = remove_mesh_face(F,v_bingren,iafacec);
%     %牙根
%     boundariestooth = select_holes_and_boundary(vertex_tooth_biaozhun,f);
%     face = fill_mesh_holes(vertex_tooth_biaozhun,f,boundariestooth,'closed',200);
%     [iaface1,ibface1]=find(bsxfun(@eq,face(:,1),index_tooth'));
%     [iaface2,ibface2]=find(bsxfun(@eq,face(iaface1,2),index_tooth'));
%     [iaface3,ibface3]=find(bsxfun(@eq,face(iaface1(iaface2),3),index_tooth'));
%     zz = 1:length(face);[iaf,ibf]=ismember(zz,iaface1(iaface2(iaface3)));
%     iaface = zz(iaf==0);
%     [face_newtooth,vertex_newtooth,fathertooth] = remove_mesh_face(face,vertex_tooth_biaozhun,iaface);
%     vertex_newtooth = (vertex_newtooth) *R + repmat((p_bingren - C),length(vertex_newtooth),1);
%     boundaries = detect_mesh_holes_and_boundary(face_new);
%     boundary = cell2mat(boundaries(1,1));
%     X_bound = V(boundary,1);
%     Y_bound = V(boundary,2);
%     Z_bound = V(boundary,3);
    
   
    
    

    
    
    

    

    

    
    
    
%     %选择适合的高度，以这个高度的点拟合平面
%     %问题的关键：高度的选择，需要是沟壑的分隔
%     %对于切牙，可以先找切缘，同样求切缘的中心点，切牙其实高度最重要
%     %问题关键，怎么找切缘
%     vertexT = up(i).vertex*up(i).orientation';
%     or = up(i).orientation*up(i).orientation';cen = up(i).center*up(i).orientation';
%     beta = 0.8;
%     zmid = min(vertexT(:,2))+beta*(max(vertexT(:,2))-min(vertexT(:,2)));%此处接收到的模型，第二列为牙齿的牙龈到咬合面的距离，即牙齿的上下位置
%     %此处需要看是上下颌牙，另，上下颌牙的不同，可能zmid值的计算中beta的选择也可能不同
%     %因为某个平面上，不一定能够有三个以上的顶点落在此平面上，所以选择定点是一个范围，这个范围的选取跟网格的精细程度有关
%     zindex = find (vertexT(:,2) >= zmid );
%     xoy = [vertexT(zindex,1),vertexT(zindex,3)];
%     modelvertexT = model2.vertices*[upax(4).x;upax(4).y;upax(4).z]';
%     
%     
%     figure()
%     plot(xoy(:,1),xoy(:,2),'r*')
%     parapoint = mean(vertexT(zindex,:));
%     
%     
%     
%     
%     %找到这个高度的中心点
%     
%     
%     %对齐，牙齿的中心点跟轴决定牙齿的对齐
%     
end


% %DBSCAN
% for i = 1:length(up)
%     
%     vertexT = [];
%     %选择适合的高度，以这个高度的点拟合平面
%     %问题的关键：高度的选择，需要是沟壑的分隔
%     %对于切牙，可以先找切缘，同样求切缘的中心点，切牙其实高度最重要
%     %问题关键，怎么找切缘
%     vertexT = up(i).vertex*up(i).orientation';
%     or = up(i).orientation*up(i).orientation';cen = up(i).center*up(i).orientation';
%     beta = 0.8;
%     zmid = min(vertexT(:,2))+beta*(max(vertexT(:,2))-min(vertexT(:,2)));%此处接收到的模型，第二列为牙齿的牙龈到咬合面的距离，即牙齿的上下位置
%     %此处需要看是上下颌牙，另，上下颌牙的不同，可能zmid值的计算中beta的选择也可能不同
%     %因为某个平面上，不一定能够有三个以上的顶点落在此平面上，所以选择定点是一个范围，这个范围的选取跟网格的精细程度有关
%     zindex = find (vertexT(:,2) >= zmid );
%     xoy = [vertexT(zindex,1),vertexT(zindex,3)];
%     modelvertexT = model2.vertices*[upax(4).x;upax(4).y;upax(4).z]';
%     
%     
% %     figure()
% %     trimesh(face,vertexT(:,1),vertexT(:,3),vertexT(:,2))
% %     hold on
% %     quiver3(cen(1),cen(3),cen(2),or(1,1),or(1,3),or(1,2),10); 
% %     quiver3(cen(1),cen(3),cen(2),or(2,1),or(2,3),or(2,2),10); 
% %     quiver3(cen(1),cen(3),cen(2),or(3,1),or(3,3),or(3,2),10); 
% %     hold off
%     figure()
%     plot(xoy(:,1),xoy(:,2),'r*')
%     parapoint = mean(vertexT(zindex,:));
%     
%     
%     
%     
%     %找到这个高度的中心点
%     
%     
%     %对齐，牙齿的中心点跟轴决定牙齿的对齐
%     
% end
% for i = 3:length(up)
%     p = [];
%     l = length(up(i).vertex);
%     face = [l+up(i).face(:,3)+1,l+up(i).face(:,2)+1,l+up(i).face(:,1)+1];
%     %提取其中的一个牙冠，建立局部坐标系，由世界坐标系转至局部坐标系
%     vertexT = up(i).vertex*up(i).orientation + up(i).center;
%     %计算平均曲率
%     [~,MC]=curvatures(vertexT(:,1),vertexT(:,3),vertexT(:,2),face);
% %     trimesh(face,vertexT(:,1),vertexT(:,3),vertexT(:,2))
% %     pp= patch('Faces',face,'Vertices',[vertexT(:,1),vertexT(:,3),vertexT(:,2)],...
% %         'FaceVertexCData',MC,'FaceColor','interp','EdgeColor','none');
% %     caxis([-1 , 0.2])
% 
% 
%     %z值替换
%     zzz = 0.5*vertexT(:,2)+0.5*MC;
%     vertexT(:,2) = zzz;
%     trimesh(face,vertexT(:,1),vertexT(:,3),vertexT(:,2))
%     %剔除低点值
%     beta = 0.625;
%     zmid = min(vertexT(:,2))+beta*(max(vertexT(:,2)) - min(vertexT(:,2)));
%     index = find(vertexT(:,2)>=zmid);
%     p = vertexT(index,:);
%     %映射（投影到XOY平面？）
%     %DBSCAN
%     epsilon=0.5;MinPts=20;
%     [IDX,~]= DBSCAN([p(:,1),p(:,3)],epsilon,MinPts);
%     PlotClusterinResult([p(:,1),p(:,3)], IDX);
%     %簇数以及中心点
%     k = max(IDX);
%     for j=1:k
%         Xi = [];
%         Xi=p(IDX==j,:);
%         center(j,:) = mean(Xi);  
%     end
%     %K-Means
%     %中心点映射
%     
%     
%     %K-Means聚类
%     
%     
%     
%     %特征点识别   
% end
