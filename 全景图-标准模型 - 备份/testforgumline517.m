clc
clear all
functionname='testforgumline517.m'; functiondir=which(functionname);
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
% % 加载标准牙
load('dental_crown.mat');
load('dentalmodelwithroot0.1forlow.mat');
load('axisofdental_crowm');
load('trfromtoothtocrown.mat');
% load('up_bingren.mat');
% load('low_bingren.mat');

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
%     addpath([functiondir 'data2\53120\u'])
%     fdi = textread('FDI__.txt');
    fdi = textread('FDIUpper__.dxt');
    n = find (fdi == x);
    if isempty(n)
        disp( '该病例没有此牙齿 ');
        return;
       
    else
        namestr1 = ['toothUpper_',num2str(n-1),'.','obj'];
        bingren = Read_Obj(namestr1);
        namestr2 = ['gumlineUpper_',num2str(n-1),'.','obj'];
        gumline =  ReadObj(namestr2);
        axis = ReadObj('AXISUpper_.obj');
        axi = axis(3*(n-1)+1:3*n,:);
   
        c = mean(bingren.vertex);y = max(gumline(:,2));
        p_bingren =[c(1),y(1),c(3)];n_bingren = axi(2,:);v_bingren = bingren.vertex;f_bingren = bingren.face;
        
        index_bingren = find(v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-2);
        point_bingren = v_bingren(index_bingren,:);
        a1 = ismember(f_bingren(:,1),index_bingren);
        aa1 = find(a1==1);
        a2 = ismember(f_bingren(:,2),index_bingren);
        aa2 = find(a2==1);
        a3 = ismember(f_bingren(:,3),index_bingren);
        aa3 = find(a3==1);
        [c1,~,~] = intersect(aa1,aa2);c2 = intersect(c1,aa3);
%         c4 = unique([aa1;aa2;aa3],'rows','stable');
%         figure()
%         trisurf(f_bingren(c2,:),v_bingren(:,1),v_bingren(:,2),v_bingren(:,3),'facecolor','c','edgecolor','b')
        %找边界
        [raw_edges_list] = query_edges_list(f_bingren(c2,:),'sorted');
        [~,~,iu] = unique(sort(raw_edges_list,2),'rows');
        [i3,i4] = histc(iu,unique(iu));
        lone_edges_idx_vect = i3(i4) == 1;
        lone_edges_list = unique(raw_edges_list(lone_edges_idx_vect,:),'rows');
%         plot_edges(v_bingren,lone_edges_list,'r','LineWidth',2);
        %边界点：
        edge_p_b_local_list = unique([lone_edges_list(:,1);lone_edges_list(:,2)]);
        %边界点在牙冠点中的位置
        [C_b,edge_p_b_local,~] = intersect(index_bingren,edge_p_b_local_list);
        edge_p_b = point_bingren(edge_p_b_local,:);
        Edge_b = zeros(size(lone_edges_list));
        for i = 1:length(edge_p_b_local)
            nu = find(lone_edges_list == edge_p_b_local_list(i));
            Edge_b(nu) = edge_p_b_local(i);
        end
        edg_b = sotr_edge(Edge_b,1);
        
              

        ff_b = f_bingren(c2,:);
        f_b0 = ones (size(ff_b));
        for j = 1:length(point_bingren)
            nu = find( ff_b == index_bingren(j));
            f_b0(sub2ind(size(f_b0), nu)) = j;
        end
        
        
%         
        
%         
%         figure()
%         plot3(point_bingren(:,1),point_bingren(:,2),point_bingren(:,3),'g*')
%         hold on

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
        index_tooth_biaozhun = find(vertex_tooth_biaozhun(:,2)<=(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) ...
                         - p_tooth_biaozhun(1))+n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - p_tooth_biaozhun(3)))/n_crown_biaozhun(2)));
        vz_tooth = abs(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(index_tooth_biaozhun,1)-p_tooth_biaozhun(1))...
                   +n_crown_biaozhun(3)*(vertex_tooth_biaozhun(index_tooth_biaozhun,3)-p_tooth_biaozhun(3)))/n_crown_biaozhun(2));
        in_tooth = find(vz_tooth == max(vz_tooth));
        point_tooth = vertex_tooth_biaozhun(index_tooth_biaozhun(in_tooth),:);
        centcrownintooth = p_crown_biaozhun+point_tooth-point_crown;%切后标准根的边缘中心
        %%利用标准牙和病人牙齿的中心点以及三轴进行对齐，对齐后按照病人牙冠的切割位置进行切割
        R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*axi;
        C=centcrownintooth*R;
        tooth_T = (vertex_tooth_biaozhun) *R + repmat((c - C),length(vertex_tooth_biaozhun),1);
        
        %整体变化牙齿的形态
        index_1 = find(v_bingren(:,2)>=(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))...
                & v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+0.5);
        point1 = v_bingren(index_1,:);
        index_2 = find(tooth_T(:,2)>=(c(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2))...
                & tooth_T(:,2)<=(p_bingren(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2))+0.5);
        point2 = tooth_T(index_2,:);
        %确定变形参数
%         point11 = point1(randperm(length(point1),length(point2)),:);
        for i = 1:length(point2)
           [minValue,r]=mindis(point1,point2(i,:),1);
           minvalue(i) = minValue;  row(i) = r;
           
        end
        point11 = point1(row,:);
%         [R,t,BRt,e,~,~] = icp(point11,point2);
%         tooth_T = bsxfun(@plus,tooth_T*R,t);
        [d,z,transform]=procrustes(point11,point2,'scaling',0,'reflection',0);
        tooth_T = transform.b * tooth_T * transform.T + transform.c(1,:);
        
        index_tooth =  find(tooth_T(:,2)>(p_bingren(2) - ( n_bingren(1)*(tooth_T(:,1) ...
                         - p_bingren(1))+ n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/ n_bingren(2))-2);%切出来的牙根的索引值%注意！按照病人牙冠方向切
        tooth = tooth_T(index_tooth ,:);
        a_t1 = ismember(f(:,1),index_tooth);
        aa_t1 = find(a_t1==1);
        a_t2 = ismember(f(:,2),index_tooth);
        aa_t2 = find(a_t2==1);
        a_t3 = ismember(f(:,3),index_tooth);
        aa_t3 = find(a_t3==1);
        [c_t1,~,~] = intersect(aa_t1,aa_t2);c_t2 = intersect(c_t1,aa_t3);
%         c_t4 = unique([aa_t1;aa_t2;aa_t3],'rows','stable');
%         figure()
%         trisurf(f_bingren(c4,:),v_bingren(:,1),v_bingren(:,2),v_bingren(:,3),'facecolor','c','edgecolor','b')
        %找边界
        [raw_edges_list_t] = query_edges_list(f(c_t2,:),'sorted');
        [~,~,iu_t] = unique(sort(raw_edges_list_t,2),'rows');
        [i3,i4] = histc(iu_t,unique(iu_t));
        lone_edges_idx_vect_t = i3(i4) == 1;
        lone_edges_list_t = unique(raw_edges_list_t(lone_edges_idx_vect_t,:),'rows');
        
        Edge_t = zeros(size(lone_edges_list_t));
%         plot_edges(tooth_T,lone_edges_list_t,'r','LineWidth',2);
       
        figure()
        trisurf(f_bingren(c2,:),v_bingren(:,1),v_bingren(:,2),v_bingren(:,3),'facecolor','c','edgecolor','b')
        hold on
        trisurf(f(c_t2,:),tooth_T(:,1),tooth_T(:,2),tooth_T(:,3),'facecolor','c','edgecolor','b')
        ff_t = f(c_t2,:);
        f_t0 = ones (size(ff_t));
        for j = 1:length(tooth)
            nu = find( ff_t == index_tooth(j));
            f_t0(sub2ind(size(f_t0), nu)) = j+length(point_bingren);
        end
        f_end = [f_b0;f_t0];v_end = [point_bingren;tooth];
%         Edge_b = sortrows(Edge_b,1);
%         figure()
%         for i = 1:length(edg_b)
%             plot_edges(v_end,edg_b(i,:),'r','LineWidth',2);
%             hold on
%         end
        

        
        %边界点：
        edge_p_t_local_list = unique([lone_edges_list_t(:,1);lone_edges_list_t(:,2)]);
        %边界点在牙根点中的位置
        [C_t,edge_p_t_local,~] = intersect(index_tooth,edge_p_t_local_list);
  
        edge_p_t = tooth(edge_p_t_local,:);
        edge_p_t_local = edge_p_t_local+length(point_bingren);
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
        tri11 = v_end(edg_b(3,1),:) - v_end(edg_b(1,1),:);tri12 = v_end(Edge_t(r,1),:) - v_end(edg_b(3,1),:);
        tri21 = v_end(Edge_t(rr,3-rm),:) - v_end(Edge_t(r,1),:);tri22 = v_end(edg_b(1,1),:) - v_end(Edge_t(rr,3-rm),:);
        s_b = sign(dot(v_end(edg_b(1,1),:) -mean(edge_p_b),cross(tri11,tri12)));
        s_t = sign(dot(v_end(Edge_t(r,1),:) - mean(edge_p_t),cross(tri21,tri22)));
        if s_b ~= s_t
           edg_t = sotr_edge(Edge_t,r);
        else
            E_t = [Edge_t(:,2),Edge_t(:,1)];
            edg_t = sotr_edge(E_t,r);
        end

        
        
        %边界点融合,以稀疏（牙根）的边线为基准，找密集（牙冠）中离两个点距离最近的点

        dt = [];
%         hold on
        for j = 1:length(edg_t)
            p1 = v_end(edg_t(j,1),:);
            p2 = v_end(edg_t(j,2),:);

            [~,ro]=mindis(v_end(edg_b(:,1),:),(p1+p2)/2,1);
            roww(j) = ro;
%             rr = edg_t(j,:);
%             rr = [find(index_tooth == lone_edges_list_t(j,1))+length(point_bingren) find(index_tooth == lone_edges_list_t(j,2))+length(point_bingren)];
%             dt = [dt;[rr,edge_p_b_local(ro)]];
            dt = [dt;[edg_t(j,:),edg_b(ro,1)]];
%             trisurf(dt(j,:),v_end(:,1),v_end(:,2),v_end(:,3),'facecolor','c','edgecolor','b')
            
        end
        if roww(1)<length(edg_b)/2
            for j = 1:length(roww)-1
                if roww(j) ~= roww(j+1) 
                    dt = [dt;[edg_b(roww(j):roww(j+1)-1,:),repmat(edg_t(j,2),roww(j+1)-roww(j),1)]];
                end
            end
            if roww(end)==length(edg_b)
                dt = [dt;[edg_b(roww(end),:),edg_t(1,1)]];
                if roww(1)~=1
                    dt = [dt;[edg_b(1:roww(1),:),repmat(edg_t(1,1),roww(1),1)]];
                end
            elseif roww(end) > length(edg_b)/2 && roww(end)<length(edg_b)
                dt = [dt;[edg_b(roww(end):length(edg_b),:),repmat(edg_t(1,1),length(edg_b)-roww(end)+1,1)]];
                if roww(1)~=1
                    dt = [dt;[edg_b(1:roww(1),:),repmat(edg_t(1,1),roww(1),1)]];
                end
            elseif  roww(end) < length(edg_b)/2 && roww(1)~=1
               s = find(roww > length(edg_b)/2);
               dt = [dt;[edg_b(roww(s(end)):length(edg_b),:),repmat(edg_t(s(end),2),length(edg_b)-roww(s(end))+1,1)]];
               dt = [dt;[edg_b(1:roww(end),:),repmat(edg_t(s(end),2),roww(end),1)]]; 
               dt = [dt;[edg_b(roww(end):roww(1)-1,:),repmat(edg_t(1,1),roww(1)-roww(end),1)]]; 
            end
            
        else
           s = find(roww<length(edg_b)/2);
%            roww = [roww roww(1:s(1)-1)];
           for j = s(1):length(roww)-1
                if roww(j) ~= roww(j+1) 
                    dt = [dt;[edg_b(roww(j):roww(j+1)-1,:),repmat(edg_t(j,2),roww(j+1)-roww(j),1)]];   
                end
                
           end 
           if s(1)>2
                for j = 1:s(1)-2
                     if roww(j) ~= roww(j+1) 
                        dt = [dt;[edg_b(roww(j):roww(j+1)-1,:),repmat(edg_t(j,2),roww(j+1)-roww(j),1)]];   
                     end
                end
           end
           if roww(end)<roww(1) 
              dt = [dt;[edg_b(roww(end):roww(1)-1,:),repmat(edg_t(end,2),roww(1)-roww(end),1)]];
           end
            if roww(s(1)-1)<=length(edg_b)
                 dt = [dt;[edg_b(roww(s(1)-1):length(edg_b),:),repmat(edg_t(s(1)-1,2),length(edg_b)-roww(s(1)-1)+1,1)]];
            end
            dt = [dt;[edg_b(1:roww(s(1))-1,:),repmat(edg_t(s(1),1),roww(s(1))-1,1)]];
            
           
        end
        ff = [f_end;dt];
        b = select_holes_and_boundary(v_end,ff);
                
%         pats = [];

        %牙冠变形控制部分
        index_crown_change = find(v_bingren(:,2)>(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-1 ...
                                & v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+2);
        crownforchange = v_bingren(index_crown_change,:);
        %牙冠牙根分别形成网格
        p = v_end;
        index_tooth_change =  find(p(:,2)>(c(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-1 ...
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
    end

else
%     addpath([functiondir 'data2\61364\l'])
%     fdi = textread('FDI__.txt');
    fdi = textread('FDILower__.dxt');
    n = find (fdi == x);
    if isempty(n)
       disp( '该病例没有此牙齿 ');
       return;
    else
        namestr1 = ['toothLower_',num2str(n-1),'.','obj'];
        bingren = Read_Obj(namestr1);
        namestr2 = ['gumlineLower_',num2str(n-1),'.','obj'];
        gumline =  ReadObj(namestr2);
        axis = ReadObj('AXISLower_.obj');
        axi = axis(3*(n-1)+1:3*n,:);
    
         %处理
        c = mean(bingren.vertex);y = min(gumline(:,2));
        p_bingren =[c(1),y(1),c(3)];n_bingren = axi(2,:);v_bingren = bingren.vertex;f_bingren = bingren.face;
        index_bingren = find(v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - ... 
            p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+0.5);
        point_bingren = v_bingren(index_bingren,:);
         a1 = ismember(f_bingren(:,1),index_bingren);
        aa1 = find(a1==1);
        a2 = ismember(f_bingren(:,2),index_bingren);
        aa2 = find(a2==1);
        a3 = ismember(f_bingren(:,3),index_bingren);
        aa3 = find(a3==1);
        [c1,~,~] = intersect(aa1,aa2);c2 = intersect(c1,aa3);
%         c4 = unique([aa1;aa2;aa3],'rows','stable');
%         figure()
%         trisurf(f_bingren(c2,:),v_bingren(:,1),v_bingren(:,2),v_bingren(:,3),'facecolor','c','edgecolor','b')
        %找边界
        [raw_edges_list] = query_edges_list(f_bingren(c2,:),'sorted');
        [~,~,iu] = unique(sort(raw_edges_list,2),'rows');
        [i3,i4] = histc(iu,unique(iu));
        lone_edges_idx_vect = i3(i4) == 1;
        lone_edges_list = unique(raw_edges_list(lone_edges_idx_vect,:),'rows');
%         plot_edges(v_bingren,lone_edges_list,'r','LineWidth',2);
        %边界点：
        edge_p_b_local_list = unique([lone_edges_list(:,1);lone_edges_list(:,2)]);
        %边界点在牙冠点中的位置
        [C_b,edge_p_b_local,~] = intersect(index_bingren,edge_p_b_local_list);
        edge_p_b = point_bingren(edge_p_b_local,:);
        Edge_b = zeros(size(lone_edges_list));
        for i = 1:length(edge_p_b_local)
            nu = find(lone_edges_list == edge_p_b_local_list(i));
            Edge_b(nu) = edge_p_b_local(i);
        end
        edg_b = sotr_edge(Edge_b,1);
        
              

        ff_b = f_bingren(c2,:);
        f_b0 = ones (size(ff_b));
        for j = 1:length(point_bingren)
            nu = find( ff_b == index_bingren(j));
            f_b0(sub2ind(size(f_b0), nu)) = j;
        end
%         figure()
%         plot3(point_bingren(:,1),point_bingren(:,2),point_bingren(:,3),'g*')
%         hold on

        %标准牙冠
        face_crown_biaozhun = crown.model.face;vertex_crown_biaozhun = crown.model.vertex;
        p_crown_biaozhun = mean(vertex_crown_biaozhun);
        n_crown_biaozhun = crown_ax.y;
        index_crown_biaozhun = find(vertex_crown_biaozhun(:,2)>=(p_crown_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_crown_biaozhun(:,1) ...
                    - p_crown_biaozhun(1))+n_crown_biaozhun(3)*(vertex_crown_biaozhun(:,3) - p_crown_biaozhun(3)))/n_crown_biaozhun(2)));
        for i =1:length(index_crown_biaozhun)
            v_crown_biaozhun(i) = abs(dot((vertex_crown_biaozhun(index_crown_biaozhun(i),:)-p_crown_biaozhun),n_crown_biaozhun)/sqrt(sum(n_crown_biaozhun.*n_crown_biaozhun)));
        end
%         v_crown_biaozhun = abs(p_crown_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_crown_biaozhun(index_crown_biaozhun,1)-p_crown_biaozhun(1))...
%                      +n_crown_biaozhun(3)*(vertex_crown_biaozhun(index_crown_biaozhun,3)-p_crown_biaozhun(3)))/n_crown_biaozhun(2));
        in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
        
        point_crown =  vertex_crown_biaozhun(index_crown_biaozhun(in_crown),:);
%         in1 = find(vertex_crown_biaozhun(:,2) == max(vertex_crown_biaozhun(:,2)));
%         point_crown =  vertex_crown_biaozhun(in1,:);

        %标准牙根
        vertex_tooth_biaozhun = root.vertices;
        f = root.faces;
        p_tooth_biaozhun = mean(vertex_tooth_biaozhun);
        index_tooth_biaozhun = find(vertex_tooth_biaozhun(:,2)>=(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) ...
                         - p_tooth_biaozhun(1))+n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - p_tooth_biaozhun(3)))/n_crown_biaozhun(2)));
        for i =1:length(index_tooth_biaozhun)
            vz_tooth(i) = abs(dot((vertex_tooth_biaozhun(index_tooth_biaozhun(i),:)-p_tooth_biaozhun),n_crown_biaozhun)/sqrt(sum(n_crown_biaozhun.*n_crown_biaozhun)));
        end
        in_tooth = find(vz_tooth == max(vz_tooth));
        
        point_tooth = vertex_tooth_biaozhun(index_tooth_biaozhun(in_tooth),:);

        vertex_tooth_biaozhun = vertex_tooth_biaozhun-repmat(point_tooth,length(vertex_tooth_biaozhun),1)+repmat(point_crown,length(vertex_tooth_biaozhun),1);
        %%利用标准牙和病人牙齿的中心点以及三轴进行对齐，对齐后按照病人牙冠的切割位置进行切割
        R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*axi;
        C=p_crown_biaozhun*R;
        tooth_T = (vertex_tooth_biaozhun) *R + repmat((c - C),length(vertex_tooth_biaozhun),1);
        figure()
        trisurf(f,tooth_T(:,1),tooth_T(:,2),tooth_T(:,3),'facecolor','c','edgecolor','b')
        hold on
        trisurf(bingren.face,v_bingren(:,1),v_bingren(:,2),v_bingren(:,3),'facecolor','c','edgecolor','r')
        
        %整体变化牙齿的形态
        index_1 = find(v_bingren(:,2)<=(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+0.5 ...
                & v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2)+0.5));
        point1 = v_bingren(index_1,:);
        index_2 = find(tooth_T(:,2)<=(c(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2))+0.5  ...
                & tooth_T(:,2)>=(p_bingren(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2)+0.5));
        point2 = tooth_T(index_2,:);
        %确定变形参数
%         point11 = point1(randperm(length(point1),length(point2)),:);
        for i = 1:length(point2)
           [minValue,r]=mindis(point1,point2(i,:),1);
           minvalue(i) = minValue;  row(i) = r;
        end
        point11 = point1(row,:);
%         [R,t,BRt,e,~,~] = icp(point11,point2);
%         tooth_T = bsxfun(@plus,tooth_T,t);
        [d,z,transform]=procrustes(point11,point2,'scaling',0,'reflection',0);
        tooth_T = bsxfun(@plus,tooth_T, transform.c(1,:));
%         tooth_T = transform.b * tooth_T * transform.T + transform.c(1,:);
        
        index_tooth =  find(tooth_T(:,2)<(p_bingren(2) - ( n_bingren(1)*(tooth_T(:,1) ...
                         - p_bingren(1))+ n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/ n_bingren(2))+0.5);%切出来的牙根的索引值%注意！按照病人牙冠方向切
        tooth = tooth_T(index_tooth ,:);
        a_t1 = ismember(f(:,1),index_tooth);
        aa_t1 = find(a_t1==1);
        a_t2 = ismember(f(:,2),index_tooth);
        aa_t2 = find(a_t2==1);
        a_t3 = ismember(f(:,3),index_tooth);
        aa_t3 = find(a_t3==1);
        [c_t1,~,~] = intersect(aa_t1,aa_t2);c_t2 = intersect(c_t1,aa_t3);
%         c_t4 = unique([aa_t1;aa_t2;aa_t3],'rows','stable');
%         figure()
%         trisurf(f_bingren(c4,:),v_bingren(:,1),v_bingren(:,2),v_bingren(:,3),'facecolor','c','edgecolor','b')
        %找边界
        [raw_edges_list_t] = query_edges_list(f(c_t2,:),'sorted');
        [~,~,iu_t] = unique(sort(raw_edges_list_t,2),'rows');
        [i3,i4] = histc(iu_t,unique(iu_t));
        lone_edges_idx_vect_t = i3(i4) == 1;
        lone_edges_list_t = unique(raw_edges_list_t(lone_edges_idx_vect_t,:),'rows');
        
        Edge_t = zeros(size(lone_edges_list_t));
        figure()
        trisurf(f_bingren(c2,:),v_bingren(:,1),v_bingren(:,2),v_bingren(:,3),'facecolor','c','edgecolor','b')
        hold on
        trisurf(f(c_t2,:),tooth_T(:,1),tooth_T(:,2),tooth_T(:,3),'facecolor','c','edgecolor','b')
        ff_t = f(c_t2,:);
        f_t0 = ones (size(ff_t));
        for j = 1:length(tooth)
            nu = find( ff_t == index_tooth(j));
            f_t0(sub2ind(size(f_t0), nu)) = j+length(point_bingren);
        end
        f_end = [f_b0;f_t0];v_end = [point_bingren;tooth];
%         Edge_b = sortrows(Edge_b,1);
%         figure()
%         for i = 1:length(edg_b)
%             plot_edges(v_end,edg_b(i,:),'r','LineWidth',2);
%             hold on
%         end
        

        
        %边界点：
        edge_p_t_local_list = unique([lone_edges_list_t(:,1);lone_edges_list_t(:,2)]);
        %边界点在牙根点中的位置
        [C_t,edge_p_t_local,~] = intersect(index_tooth,edge_p_t_local_list);
  
        edge_p_t = tooth(edge_p_t_local,:);
        edge_p_t_local = edge_p_t_local+length(point_bingren);
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
        tri11 = v_end(edg_b(3,1),:) - v_end(edg_b(1,1),:);tri12 = v_end(Edge_t(r,1),:) - v_end(edg_b(3,1),:);
        tri21 = v_end(Edge_t(rr,3-rm),:) - v_end(Edge_t(r,1),:);tri22 = v_end(edg_b(1,1),:) - v_end(Edge_t(rr,3-rm),:);
        s_b = sign(dot(v_end(edg_b(1,1),:) -mean(edge_p_b),cross(tri11,tri12)));
        s_t = sign(dot(v_end(Edge_t(r,1),:) - mean(edge_p_t),cross(tri21,tri22)));
        if s_b ~= s_t
           edg_t = sotr_edge(Edge_t,r);
        else
            E_t = [Edge_t(:,2),Edge_t(:,1)];
            edg_t = sotr_edge(E_t,r);
        end

        
        
        %边界点融合,以稀疏（牙根）的边线为基准，找密集（牙冠）中离两个点距离最近的点

        dt = [];
%         hold on
        for j = 1:length(edg_t)
            p1 = v_end(edg_t(j,1),:);
            p2 = v_end(edg_t(j,2),:);

            [~,ro]=mindis(v_end(edg_b(:,1),:),(p1+p2)/2,1);
            roww(j) = ro;
%             rr = edg_t(j,:);
%             rr = [find(index_tooth == lone_edges_list_t(j,1))+length(point_bingren) find(index_tooth == lone_edges_list_t(j,2))+length(point_bingren)];
%             dt = [dt;[rr,edge_p_b_local(ro)]];
            dt = [dt;[edg_t(j,:),edg_b(ro,1)]];
%             trisurf(dt(j,:),v_end(:,1),v_end(:,2),v_end(:,3),'facecolor','c','edgecolor','b')
            
        end
        if roww(1)<length(edg_b)/2
            for j = 1:length(roww)-1
                if roww(j) ~= roww(j+1) 
                    dt = [dt;[edg_b(roww(j):roww(j+1)-1,:),repmat(edg_t(j,2),roww(j+1)-roww(j),1)]];
                end
            end
            if roww(end)==length(edg_b)
                dt = [dt;[edg_b(roww(end),:),edg_t(1,1)]];
                if roww(1)~=1
                    dt = [dt;[edg_b(1:roww(1),:),repmat(edg_t(1,1),roww(1),1)]];
                end
            elseif roww(end) > length(edg_b)/2 && roww(end)<length(edg_b)
                dt = [dt;[edg_b(roww(end):length(edg_b),:),repmat(edg_t(1,1),length(edg_b)-roww(end)+1,1)]];
                if roww(1)~=1
                    dt = [dt;[edg_b(1:roww(1),:),repmat(edg_t(1,1),roww(1),1)]];
                end
            elseif  roww(end) < length(edg_b)/2 && roww(1)~=1
               s = find(roww > length(edg_b)/2);
               dt = [dt;[edg_b(roww(s(end)):length(edg_b),:),repmat(edg_t(s(end),2),length(edg_b)-roww(s(end))+1,1)]];
               dt = [dt;[edg_b(1:roww(end),:),repmat(edg_t(s(end),2),roww(end),1)]]; 
               dt = [dt;[edg_b(roww(end):roww(1)-1,:),repmat(edg_t(1,1),roww(1)-roww(end),1)]]; 
            end
            
        else
           s = find(roww<length(edg_b)/2);
%            roww = [roww roww(1:s(1)-1)];
           for j = s(1):length(roww)-1
                if roww(j) ~= roww(j+1) 
                    dt = [dt;[edg_b(roww(j):roww(j+1)-1,:),repmat(edg_t(j,2),roww(j+1)-roww(j),1)]];   
                end
                
           end 
           if s(1)>2
                for j = 1:s(1)-2
                     if roww(j) ~= roww(j+1) 
                        dt = [dt;[edg_b(roww(j):roww(j+1)-1,:),repmat(edg_t(j,2),roww(j+1)-roww(j),1)]];   
                     end
                end
           end
           if roww(end)<roww(1) 
              dt = [dt;[edg_b(roww(end):roww(1)-1,:),repmat(edg_t(end,2),roww(1)-roww(end),1)]];
           end
            if roww(s(1)-1)<=length(edg_b)
                 dt = [dt;[edg_b(roww(s(1)-1):length(edg_b),:),repmat(edg_t(s(1)-1,2),length(edg_b)-roww(s(1)-1)+1,1)]];
            end
            dt = [dt;[edg_b(1:roww(s(1))-1,:),repmat(edg_t(s(1),1),roww(s(1))-1,1)]];
            
           
        end
        ff = [f_end;dt];
        p = v_end;
        b = select_holes_and_boundary(v_end,ff);
        
        
        
        
%         plot3(tooth(:,1),tooth(:,2),tooth(:,3),'r*')

%         tooth = vertex_tooth_biaozhun(index_tooth,:);
%         R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*bingren.orientation;
%         C=centcrownintooth*R;
%         tooth_T = (tooth) *R + repmat((p_bingren - C),length(tooth),1);

        %牙冠变形控制部分
        index_crown_change = find(v_bingren(:,2)<(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+1 ...
                                & v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-2);
        crownforchange = v_bingren(index_crown_change,:);
        
        
        
        

        index_tooth_change =  find(p(:,2)<(c(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+1 ...
            &p(:,2)>=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-2);
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
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+1.5 ...
            |p(:,2)<=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-3);
    end
end
    
    
    


[Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(p_T,1), ff, exterior);
[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(p,ff, bi_bndtype,masstype,reduction,Omega,N0,N1);
%找到牙根部分跟牙冠变形部分最近的点，将牙根上的点拉至牙冠处
for i = 1:length(pchange)
   [minValue,r]=mindis(crownforchange,pchange(i,:),1);
   minvalue(i) = minValue;  row(i) = r;
   p_T(index_tooth_change(i),:) = crownforchange(r,:);
end
bi_V = biharm_solve_with_factor( ...
    bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
    ff, p_T, Omega, N0, N1, bi_bndtype, reduction,BZ1,p);
figure()
trisurf(ff,bi_V(:,1),bi_V(:,2),bi_V(:,3),'facecolor','c','edgecolor','b')
% axis image
namestr3 = [num2str(x),'.','obj'];
writeOBJ(namestr3,bi_V,ff)


 

writeOBJ('最终网格.obj',bi_V, ff)
disp('finished.');
 
