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
    %     if s == 2
    %        bingren = up(8+g); 
    %     end
    %     if s == 1
    %        bingren = up(9-g); 
    %     end
        %处理
        c = mean(bingren.vertex);y = max(gumline(:,2));
        p_bingren =[c(1),y(1),c(3)];n_bingren = axi(2,:);v_bingren = bingren.vertex;f_bingren = bingren.face;
%         b = select_holes_and_boundary(v_bingren,f_bingren);
%         ff = fill_mesh_holes(v_bingren,f_bingren,b,'closed',200);
        
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
        [R,t,BRt,e,~,~] = icp(point11,point2);
        tooth_T = bsxfun(@plus,tooth_T*R,t);
%         [d,z,transform]=procrustes(point11,point2,'scaling',0,'reflection',0);
%         tooth_T = transform.b * tooth_T * transform.T + transform.c(1,:);
        
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
%         q1 = v_end(edg_b(1,1),:) - v_end(edg_b(1,2),:);q2 = v_end(Edge_t(r,1),:) - v_end(Edge_t(r,2),:);
%         s_b = q1 - dot(axi(2,:),q1)/norm(axi(2,:), 2)*axi(2,:);
%         s_t = q2 - dot(axi(2,:),q2)/norm(axi(2,:), 2)*axi(2,:);
%         s_b = sign(dot(axi(2,:),cross(v_end(edg_b(1,1),:) -mean(edge_p_b) , v_end(edg_b(1,2),:) - mean(edge_p_b))));
%         s_t = sign(dot(axi(2,:),cross(v_end(Edge_t(r,1),:) -mean(edge_p_t) , v_end(Edge_t(r,2),:) - mean(edge_p_t))));
        if s_b ~= s_t
           edg_t = sotr_edge(Edge_t,r);
        else
            E_t = [Edge_t(:,2),Edge_t(:,1)];
            edg_t = sotr_edge(E_t,r);
        end
%         p1 = v_end(unique([edg_t(:,1);edg_t(:,2)]),:);
%         sp = spline(1:length(edg_t),[p1(:,1)';p1(:,2)';p1(:,3)'],1:length(edg_b))';
        
%         f_v = [];
%         n = round(length(edg_b)/length(edg_t));
%         if length(edg_b)/n > length(edg_t)
%             for i = 1:length(edg_t)-1
%                 f_v = [f_v;[edg_b((i-1)*n+1:i*n,:),repmat(edg_t(i,1),n,1)]];
%                 f_v = [f_v;[edg_t(i,:),edg_b(i*n,2)]];
%                 t = i*n+1;
%             end
%             f_v = [f_v;[edg_b(t:end,:),repmat(edg_t(end,1),size(edg_b(t:end,:),1),1)]];
%             f_v = [f_v;[edg_t(end,:),edg_b(1,1)]];
%         elseif length(edg_b)/n == length(edg_t)
%             for i = 1:length(edg_t)-1
%                 f_v = [f_v;[edg_b((i-1)*n+1:i*n,:),repmat(edg_t(i,1),n,1)]];
%                 f_v = [f_v;[edg_t(i,:),edg_b(i*n,2)]];
%                 t = i*n+1;
%             end
%             f_v = [f_v;[edg_t(end,:),edg_b(1,1)]];
%         else
%             for i = 1:length(edg_b)/n-1
%                 f_v = [f_v;[edg_b((i-1)*n+1:i*n,:),repmat(edg_t(i,1),n,1)]];
%                 f_v = [f_v;[edg_t(i,:),edg_b(i*n,2)]];
%                 t = i*n+1;
%             end
%             f_v = [f_v;[edg_t(i+1:end,:),repmat(edg_b(end,2),size(edg_t(i+1:end,:),1),1)]];
%             f_v = [f_v;[edg_b(t:end,:),repmat(edg_t(end,1),size(edg_b(t:end,:),1),1)]];
%         end
%         ff = [f_end;f_v];
            
        
%         
%         
%         edg_point_b = v_end(Edge_b(:,1),:);
%        
%         for i = 1:length(Edge_t)
%             p_t = v_end(Edge_t(i,1),:);
%             [~,ro]=mindis(edg_point_b,p_t,1);
%             ros(i) = ro;
%         end
%         EE = sortrows([Edge_t,ros'],3);
%         for i = 1:length(Edge_t)-1
%             if EE(i+1,3) == EE(i,3)
%                 fv = [f_v;[EE(i,1:2),Edge_b(EE(i,3),1)];[EE(i+1,1:2),Edge_b(EE(i+1,3),1)]];
%             else
%                 f_v = [f_v;[EE(i,1:2),Edge_b(EE(i+1,3),1)]];
%                 f_v = [f_v;[Edge_b(EE(i,3):EE(i+1,3),:),repmat(EE(i,1),size(Edge_b(EE(i,3):EE(i+1,3),:),1),1)]];%以牙冠上的线为边，牙根上的点
%             end
%         end
%         
        
        
        
        %边界点融合,以稀疏（牙根）的边线为基准，找密集（牙冠）中离两个点距离最近的点

        dt = [];
%         hold on
        for j = 1:length(edg_t)
            p1 = v_end(edg_t(j,1),:);
            p2 = v_end(edg_t(j,2),:);
%             p1 = tooth(find(index_tooth == lone_edges_list_t(j,1)),:);
%             p2 = tooth(find(index_tooth == lone_edges_list_t(j,2)),:);
            
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
            if roww(end) > length(edg_b)/2 && roww(end)<length(edg_b)
                dt = [dt;[edg_b(roww(end):length(edg_b),:),repmat(edg_t(1,1),length(edg_b)-roww(end),1)]];
            elseif  roww(end) < length(edg_b)/2
               dt = [dt;[edg_b(1:roww(end),:),repmat(edg_t(end,1),roww(end),1)]]; 
            end
            if roww(1)~=1
                dt = [dt;[edg_b(1:roww(1)-1,:),repmat(edg_t(1,1),roww(1)-1,1)]];
            end
        else
           s = find(roww<length(edg_b)/2);
           roww = [roww roww(1:s(1)-1)];
           for j = s(1):length(roww)-1
                if roww(j) ~= roww(j+1) 
                    dt = [dt;[edg_b(roww(j):roww(j+1),:),repmat(edg_t(j,2),roww(j+1)-roww(j)+1,1)]];   
                end
                
           end 
           dt = [dt;[edg_b(1:roww(s(1))-1,:),repmat(edg_t(1,2),roww(s(1))-1,1)]];
%            if r
%                dt = [dt;[edg_b(1:s(1)-1,:),repmat(edg_t(1,2),s(1)-1,1)]];
%            end
        end
        b = select_holes_and_boundary(v_end,[f_end;dt]);
                
        pats = [];
%         for j = 1:length(dt)
%             [nn,mm] = ind2sub(size(dt(:,1:2)),find(dt(:,1:2) == dt(j,1)));
%             for i = 1:length(nn)
%                 if nn(i) ~=j
%                    if dt(j,3)~=dt(nn(i),3)
%                        pat = [dt(j,1),dt(j,3),dt(nn(i),3)];
%                    else
%                        pat = [];
%                    end
%                    pats = [pats;pat];
%                 end
%             end
%         end
%         paa = unique([dt;pats],'rows');
%         f_t = [f_end;paa];
        %5.19日
        %1.将牙冠跟牙根的点分别排序，牙根点的起始点为与牙冠起始点最近的点，确保牙根与牙冠的方向一致
        %牙冠中检测到边lone_edges_list跟paa中不重合的,不重合的就跟
%         EE = unique([dt(:,1:2);dt(:,2:3);[dt(:,1) dt(:,3)]],'rows');
%         [~,~,ib] = intersect(EE,lone_edges_list,'rows');
        
        
        
        



        %更新边界
        
        
%         f_v = [];
%         [p,ia,ib] = intersect(dt(:,3),edg_b(:,1));
%         mix = sortrows([p,ia,ib],3);
%         for j = 1:length(p)-1
%             [c,~,~] = intersect(dt(mix(j,2),1:2)',dt(mix(j+1,2),1:2));
%             if length(c) == 1
%                 lines = edg_b(mix(j,3):mix(j+1,3),:);
%                 f_v = [f_v;[lines,repmat(c,size(lines,1),1)]];
%             else
%                 n1 = find(dt(:,3) == dt(mix(j,2),3));
%                 n2 = find(dt(:,3) == dt(mix(j+1,2),3));
%                 if ~isempty(n1)
%                     [c,~,~] = intersect(dt(n1,1:2),dt(mix(j+1,2),1:2));
%                      if length(c) == 1
%                         lines = edg_b(mix(j,3):mix(j+1,3),:);
%                         f_v = [f_v;[lines,repmat(c,size(lines,1),1)]];
%                      end
%                 elseif ~isempty(n2)
%                     [c,~,~] = intersect(dt(n2,1:2),dt(mix(j+1,2),1:2));
%                     if length(c) == 1
%                         lines = edg_b(mix(j,3):mix(j+1,3),:);
%                         f_v = [f_v;[lines,repmat(c,size(lines,1),1)]];
%                     end
%                 end
% 
%             end
% 
%         end
%         
%         f1 = [dt;f_v];
%             
%             if length(p) == 2
%                 [c,~,~] = intersect(dt(ia(1),1:2),dt(ia(2),1:2));
%                 if c == 1
%                    f_v = [f_v;[vp1',c]];
% %                 else
% %                     p1 = point_bingren(vp1(1),:);
% %                     p2 = point_bingren(vp1(2),:);
% %                     [~,ro]=mindis(edge_p_t,(p1+p2)/2,1);
% %                     
% %                     f_v = [f_v;[vp1',edge_p_t_local(ro)]];
%                 end
%             elseif length(p) == 1
%                 %看dt在牙根上的两个点哪个离另一个更近
%                 point = point_bingren(vp1(find([1,2]~=ib)),:);
%                 [~,ro]=mindis(v_end(dt(ia,1:2)',:),point,1);
%                 f_v = [f_v;[dt(ia,ro),vp1']];
%             else
%                 p1 = point_bingren(vp1(1),:);
%                 p2 = point_bingren(vp1(2),:);
%                 [~,ro1]=mindis(v_end(dt(:,3),:),p1,1);
%                 [~,ro2]=mindis(v_end(dt(:,3),:),p2,1);
%                 if ro1 == ro2
%                     [~,ro3]=mindis(v_end(dt(ro1,1:2)',:),p2,1);
%                     f_v = [f_v;[vp1',dt(ro1,ro3)]];
%                 else
%                    [c,~,~] = intersect((dt(ro1,1:2)),dt(ro2,1:2));
%                    if isempty(c) == 0
%                        f_v = [f_v;[vp1',c]];
% %                    else
% %                        [~,ro]=mindis(edge_p_t,(p1+p2)/2,1);
% %                        f_v = [f_v;[vp1',edge_p_t_local(ro)]];
%                    end
%                 end
%                     
%             end
%         end
                
            
        
        
        
%         %5.18思路
%         [holeCellArray,bounding_triangles,~] = findTriMeshHoles(f_t,v_end);
%         for i = 1:length(holeCellArray)
%             arr = unique(holeCellArray{i})';
%             if length(arr) == 2
%                 [~,ro_t]=mindis(tooth,mean(arr),3);
%                 lo = find(ismember(ro_t) == 0);
%                 f_t = [f_t;[arr,ro_t(lo(1))]];
%             elseif length(arr) == 3
%                 f_t = [f_t;arr];
%             else
%                 for j = 1:length(arr)-2
%                     f_t = [f_t;arr(j:j+2)];
%                 end
%             end
%         end
%         
%         
%         
%         for i = 1:length(holeCellArray)
%             arr = unique(holeCellArray{i})';
%             if length(arr) == 2
%                 [~,ro_t]=mindis(tooth,mean(arr),3);
%                 lo = find(ismember(ro_t) == 0);
%                 f_t = [f_t;[arr,ro_t(lo(1))]];
%             elseif length(arr) == 3
%                 f_t = [f_t;arr];
%             elseif length(arr) == 4
%                 f_t = [f_t;arr(1:3);arr(2:4)];
%             else 
%                 nn = find(arr>length(point_bingren));
%                 mm = find(arr<=length(point_bingren));
%                 if length(nn)>=2
%                     if mm == 0
%                       
%                       for j = 1:length(arr)-2
%                           f_t = [f_t;arr(i:i+2)];
%                       end
%                     elseif mm ==1
%                         ar1 = arr(nn);
%                         for j = 1:length(ar1)-1
%                           f_t = [f_t;[ar1(i:i+1),arr(mm)]]
%                         end
%                     else
%                        for i = 1:length(mm)-1
%                            pos = [v_end(arr(mm(i)),:);v_end(arr(mm(i+1)),:)];
%                            [~,ro_t]=mindis(v_end(arr(nn),:),mean(pos),1);
%                            f_t = [f_t;[arr(mm(i)),arr(mm(i+1)),arr(nn(ro_t))]];
%                        end
%                     end
%                 elseif length(nn)== 1
%                     ar1 = arr(mm);
%                     for j = 1:length(ar1)-1
%                        f_t = [f_t;[ar1(j:j+1),arr(nn)]];
%                     end
%                 else
%                    for j = 1:length(arr)-2
%                        f_t = [f_t;arr(j:j+2)];
%                    end
%                    
%                 end
%  
%             end
% 
%         end
%             
%         b = select_holes_and_boundary(v_end,f_t);
%         for k = 1:length(b)
%             if length(b{k}) == 3
%                 f_t = [f_t;b{k}];
%             end
%         end
%         
%         ttt = [f_end;dt];
%         %边界点融合,以密集（牙冠）的边线为基准，找稀疏（牙根）中离两个点距离最近的点
%         dtt =[];
%          for j = 1:length(lone_edges_list)
%             p1 = point_bingren(find(index_bingren == lone_edges_list(j,1)),:);
%             p2 = point_bingren(find(index_bingren == lone_edges_list(j,2)),:);    
%             [~,ro]=mindis(edge_p_t,(p1+p2)/2,1);
%             roww(j) = ro;
%             rr = [find(index_bingren == lone_edges_list(j,1)) find(index_bingren == lone_edges_list(j,2))];
%             dtt = [dtt;[rr,edge_p_t_local(ro)]];
%             
%          end
%            patts = [];
%         for j = 1:length(dtt)
%             [nn,mm] = ind2sub(size(dtt(:,1:2)),find(dtt(:,1:2) == dtt(j,1)));
%             for i = 1:length(nn)
%                 if nn(i) ~=j
%                    if dtt(j,3)~=dtt(nn(i),3)
%                        patt = [dtt(j,1),dtt(j,3),dtt(nn(i),3)];
%                    else
%                        patt = [];
%                    end
%                    patts = [patts;patt];
%                 end
%             end
%         end
%         paat = unique([dtt;patts],'rows');
%         f_tt = [f_end;paat];
%         b = select_holes_and_boundary(v_end,f_tt);
%         for k = 1:length(b)
%             if length(b{k}) == 3
%                 f_tt = [f_tt;b{k}];
%             end
%         end
        
         
%         figure()
%         trisurf(dtt,v_end(:,1),v_end(:,2),v_end(:,3),'facecolor','c','edgecolor','b')
%         tt = [f_end;dtt];
% %         lines = [];patchs = [];
% %         for i = 1:length(edge_p_t)
% %             [~,r1]=mindis(edge_p_b,edge_p_t(i,:),1);
% %             line1 = [edge_p_t_local(i),edge_p_b_local(r1)];
% %             [~,r2]=mindis(edge_p_b,edge_p_b(r1,:),2);
% %             patch1 = [edge_p_t_local(i),edge_p_b_local(r1),edge_p_b_local(r2(2))];
% %             [~,r3]=mindis(edge_p_t,edge_p_b(r2(2),:),2);
% %             if r3(1) ~=i
% %                 patch2 = [edge_p_t_local(i),edge_p_b_local(r2(2)),edge_p_t_local(r3(1))];
% %             else
% %                 patch2 = [edge_p_t_local(i),edge_p_b_local(r2(2)),edge_p_t_local(r3(2))];
% %             end
% %             lines = [lines;line1];patchs = [patchs;patch1;patch2];
%         end
        
       
        
        
        
%         %未完待续5.17
%         for i = 1:length(roww)
%             [k1,k2] = ind2sub(size(lone_edges_list),find(lone_edges_list == roww));
%             for j =1:length(k1)
%                 linee = lone_edges_list(k1(j),:);
%                 p1 = point_bingren(find(index_bingren == lone_edges_list(k1(j))),:);
%                 p2 = point_bingren(find(index_bingren == lone_edges_list(k1(j))),:);
%                 pp1 = tooth(find(index_tooth == lone_edges_list_t(j,1)),:);
%                 pp2 = tooth(find(index_tooth == lone_edges_list_t(j,2)),:);
%             end
%                 
%             
%         end
            
        
            
       
        
        
        
        
        
        
        
  
%         plot3(tooth(:,1),tooth(:,2),tooth(:,3),'r*')

%         tooth = vertex_tooth_biaozhun(index_tooth,:);
%         R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*bingren.orientation;
%         C=centcrownintooth*R;
%         tooth_T = (tooth) *R + repmat((p_bingren - C),length(tooth),1);

        %牙冠变形控制部分
        index_crown_change = find(v_bingren(:,2)>(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-1 ...
                                & v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+2);
        crownforchange = v_bingren(index_crown_change,:);
        %牙冠牙根分别形成网格
%         
%         [tooth_t]=MyCrustOpen(point_bingren);
%         
%         figure()
%         trimesh(tooth_t,point_bingren(:,1),point_bingren(:,2),point_bingren(:,3))
%         hold on
%         [tooth_t2]=MyCrustOpen(tooth); 
%         trimesh(tooth_t2,tooth(:,1),tooth(:,2),tooth(:,3))
%         hold off
%         %牙冠和牙根形成网格
%         p = [point_bingren;tooth];
%         [t]=MyCrustOpen(p);
%         figure()
%         trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
%         %网格补洞
%         b = select_holes_and_boundary(p,t);
%         ff = fill_mesh_holes(p,t,b,'closed',200);
%         ff = double(ff);
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
        p_bingren =[c(1),y(1),c(3)];n_bingren = axi(2,:);v_bingren = bingren.vertex;
        index_bingren = find(v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+1);
        point_bingren = v_bingren(index_bingren,:);
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
%         vz_tooth = abs(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(index_tooth_biaozhun,1)-p_tooth_biaozhun(1))...
%                    +n_crown_biaozhun(3)*(vertex_tooth_biaozhun(index_tooth_biaozhun,3)-p_tooth_biaozhun(3)))/n_crown_biaozhun(2));
        for i =1:length(index_tooth_biaozhun)
            vz_tooth(i) = abs(dot((vertex_tooth_biaozhun(index_tooth_biaozhun(i),:)-p_tooth_biaozhun),n_crown_biaozhun)/sqrt(sum(n_crown_biaozhun.*n_crown_biaozhun)));
        end
        in_tooth = find(vz_tooth == max(vz_tooth));
        
        point_tooth = vertex_tooth_biaozhun(index_tooth_biaozhun(in_tooth),:);
%         in2 = find(vertex_tooth_biaozhun(:,2) == max(vertex_tooth_biaozhun(:,2)));
%         point_tooth =  vertex_tooth_biaozhun(in2,:);
%         centcrownintooth = p_crown_biaozhun+point_tooth-point_crown;%切后标准根的边缘中心
        vertex_tooth_biaozhun = vertex_tooth_biaozhun-point_tooth+point_crown;
        %%利用标准牙和病人牙齿的中心点以及三轴进行对齐，对齐后按照病人牙冠的切割位置进行切割
        R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*axi;
        C=p_crown_biaozhun*R;
        tooth_T = (vertex_tooth_biaozhun) *R + repmat((c - C),length(vertex_tooth_biaozhun),1);
        figure()
        trisurf(f,tooth_T(:,1),tooth_T(:,2),tooth_T(:,3),'facecolor','c','edgecolor','b')
        hold on
        trisurf(bingren.face,v_bingren(:,1),v_bingren(:,2),v_bingren(:,3),'facecolor','c','edgecolor','r')
        
        %整体变化牙齿的形态
        index_1 = find(v_bingren(:,2)<=(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2)) ...
                & v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+0.5);
        point1 = v_bingren(index_1,:);
        index_2 = find(tooth_T(:,2)<=(c(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2)) ...
                & tooth_T(:,2)>=(p_bingren(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2))+0.5);
        point2 = tooth_T(index_2,:);
        %确定变形参数
%         point11 = point1(randperm(length(point1),length(point2)),:);
        for i = 1:length(point2)
           [minValue,r]=mindis(point1,point2(i,:),1);
           minvalue(i) = minValue;  row(i) = r;
        end
        point11 = point1(row,:);
        [R,t,BRt,e,~,~] = icp(point11,point2);
        tooth_T = bsxfun(@plus,tooth_T,t);
%         [d,z,transform]=procrustes(point11,point2,'scaling',0,'reflection',0);
%         tooth_T = transform.b * tooth_T * transform.T + transform.c(1,:);
        
        index_tooth =  find(tooth_T(:,2)<(p_bingren(2) - ( n_bingren(1)*(tooth_T(:,1) ...
                         - p_bingren(1))+ n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/ n_bingren(2))+1);%切出来的牙根的索引值%注意！按照病人牙冠方向切
        tooth = tooth_T(index_tooth ,:);
        
        
        
        
        
        
%         plot3(tooth(:,1),tooth(:,2),tooth(:,3),'r*')

%         tooth = vertex_tooth_biaozhun(index_tooth,:);
%         R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*bingren.orientation;
%         C=centcrownintooth*R;
%         tooth_T = (tooth) *R + repmat((p_bingren - C),length(tooth),1);

        %牙冠变形控制部分
        index_crown_change = find(v_bingren(:,2)<(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+1 ...
                                & v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-1);
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
        p = [point_bingren;tooth];
        [t]=MyCrustOpen(p);
        figure()
        trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
        %网格补洞
        b = select_holes_and_boundary(p,t);
        ff = fill_mesh_holes(p,t,b,'closed',200);
        ff = double(ff);
        index_tooth_change =  find(p(:,2)<(c(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+1 ...
            &p(:,2)>=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-1);
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
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+1 ...
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

