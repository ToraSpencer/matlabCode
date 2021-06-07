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
addpath([functiondir '������루��׼���Ͳ������ڶ��룩'])
addpath([functiondir 'siofmodeltooth'])
addpath([functiondir 'attachments'])
addpath([functiondir 'modelStadard'])
addpath([functiondir 'testdata'])
addpath([functiondir 'mesh process'])
addpath([functiondir 'testdata(4)'])
% % ���ر�׼��
load('dental_crown.mat');
load('dentalmodelwithroot0.1forlow.mat');
load('axisofdental_crowm');
load('trfromtoothtocrown.mat');
 
 

 
x = 41;
for i = 1:28
    if  dentalwithtooth1(i).ID == x
        rootTooth = dentalwithtooth1(i);
    end
end
 
for j = 1:28
    if (upax(j).fid == x)
        axisStandard = upax(j);
    end
end
for k = 1:28
    if (dental_crown(k).fid == x)
        crown = dental_crown(k);
    end
end
crownTooth = crown.model;          % ��׼����
s = fix(x/10);  %ȡ��
g = mod(x,10);%ȡ��
if s ==1 || s ==2
 
else
    fdi = textread('FDILower__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
       disp( '�ò���û�д����� ');
       return;
    else
        
        %%  1.�и������
        namestr1 = ['toothLower_',num2str(toothIdx-1),'.','obj'];
        patientTooth = Read_Obj(namestr1);
        namestr2 = ['gumlineLower_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);
        axis = ReadObj('AXISLower_.obj');
        axisPatient = axis(3*(toothIdx-1)+1:3*toothIdx,:);
    
         %����
        centerPatient = mean(patientTooth.vertex);y = min(gumline(:,2));
        p_patient =[centerPatient(1),y(1),centerPatient(3)];patientYdir = axisPatient(2,:);patientTooth.vertex = patientTooth.vertex;patientTooth.face = patientTooth.face;
        cutPatientIdx = find(patientTooth.vertex(:,2)>=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - ... 
            p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+0.5);
        cutPatientCrownVers = patientTooth.vertex(cutPatientIdx,:);
         a1 = ismember(patientTooth.face(:,1),cutPatientIdx);
        aa1 = find(a1==1);
        a2 = ismember(patientTooth.face(:,2),cutPatientIdx);
        aa2 = find(a2==1);
        a3 = ismember(patientTooth.face(:,3),cutPatientIdx);
        aa3 = find(a3==1);
        [c1,~,~] = intersect(aa1,aa2);triIdxInCutPatient = intersect(c1,aa3);

        %�ұ߽�
        [raw_edges_list] = query_edges_list(patientTooth.face(triIdxInCutPatient,:),'sorted');
        [~,~,iu] = unique(sort(raw_edges_list,2),'rows');
        [i3,i4] = histc(iu,unique(iu));
        lone_edges_idx_vect = i3(i4) == 1;
        lone_edges_list = unique(raw_edges_list(lone_edges_idx_vect,:),'rows');
        
 
        %�߽�㣺
        edge_p_b_local_list = unique([lone_edges_list(:,1);lone_edges_list(:,2)]);
        %�߽�������ڵ��е�λ��
        [C_b,edge_p_b_local,~] = intersect(cutPatientIdx,edge_p_b_local_list);
        edge_p_b = cutPatientCrownVers(edge_p_b_local,:);
        Edge_b = zeros(size(lone_edges_list));
        for i = 1:length(edge_p_b_local)
            nu = find(lone_edges_list == edge_p_b_local_list(i));
            Edge_b(nu) = edge_p_b_local(i);
        end
        patientEdge = sotr_edge(Edge_b,1);
        
             
        ff_b = patientTooth.face(triIdxInCutPatient,:);
        f_b0 = ones (size(ff_b));
        for j = 1:length(cutPatientCrownVers)
            nu = find( ff_b == cutPatientIdx(j));
            f_b0(sub2ind(size(f_b0), nu)) = j;
        end

        
        
%%
% 2. ���롪�����ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж���


%  2.1 ȷ����һ����תƽ���������     
        centerCrown = mean(crownTooth.vertex);
        rootYdir = axisStandard.y;
        index_crown_biaozhun = find(crownTooth.vertex(:,2)>=(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(:,1) ...
                    - centerCrown(1))+rootYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/rootYdir(2)));
       
        for i =1:length(index_crown_biaozhun)
            v_crown_biaozhun(i) = abs(dot((crownTooth.vertex(index_crown_biaozhun(i),:)-centerCrown),rootYdir)/sqrt(sum(rootYdir.*rootYdir)));
        end
        
        in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
        point_crown =  crownTooth.vertex(index_crown_biaozhun(in_crown),:);


        %��׼����
        p_tooth_biaozhun = mean(rootTooth.vertices);
        index_tooth_biaozhun = find(rootTooth.vertices(:,2)>=(p_tooth_biaozhun(2) - (rootYdir(1)*(rootTooth.vertices(:,1) ...
                         - p_tooth_biaozhun(1))+rootYdir(3)*(rootTooth.vertices(:,3) - p_tooth_biaozhun(3)))/rootYdir(2)));
        
        for i =1:length(index_tooth_biaozhun)
            vz_tooth(i) = abs(dot((rootTooth.vertices(index_tooth_biaozhun(i),:)-p_tooth_biaozhun),rootYdir)/sqrt(sum(rootYdir.*rootYdir)));
        end
        in_tooth = find(vz_tooth == max(vz_tooth));
        
        point_tooth = rootTooth.vertices(index_tooth_biaozhun(in_tooth),:);

        rootTooth.vertices = rootTooth.vertices-repmat(point_tooth,length(rootTooth.vertices),1)+repmat(point_crown,length(rootTooth.vertices),1);
        
        
%       2.2 ��һ����תƽ��
temp = inv([axisStandard.x;axisStandard.y;axisStandard.z]);
R = temp * axisPatient;
newCent = centerCrown*R;
        movedRootTooth2 = (rootTooth.vertices) *R + repmat((centerPatient - newCent),length(rootTooth.vertices),1);
 
        %����仯���ݵ���̬
        index_1 = find(patientTooth.vertex(:,2)<=(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2)) ...
                & patientTooth.vertex(:,2)>=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2)+0.5));
        point1 = patientTooth.vertex(index_1,:);
        index_2 = find(movedRootTooth2(:,2)<=(centerPatient(2) - (patientYdir(1)*(movedRootTooth2(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth2(:,3) - p_patient(3)))/patientYdir(2))  ...
                & movedRootTooth2(:,2)>=(p_patient(2) - (patientYdir(1)*(movedRootTooth2(:,1) - p_patient(1))+patientYdir(3)*(movedRootTooth2(:,3) - p_patient(3)))/patientYdir(2)+0.5));
        point2 = movedRootTooth2(index_2,:);
        %ȷ�����β���

        for i = 1:length(point2)
           [minValue,startEdgeIdx]=mindis(point1,point2(i,:),1);
           minvalue(i) = minValue;  row(i) = startEdgeIdx;
        end
        point11 = point1(row,:);
        [R,t,BRt,e,~,~] = icp(point11,point2);
        movedRootTooth2 = bsxfun(@plus,movedRootTooth2,t);

        
%% 3. �и������׼��  
%       3.1 ȷ���и�㡣
        rootCutIdx =  find(movedRootTooth2(:,2)<(p_patient(2) - ( patientYdir(1)*(movedRootTooth2(:,1) ...
                         - p_patient(1))+ patientYdir(3)*(movedRootTooth2(:,3) - p_patient(3)))/ patientYdir(2))+0.5);%�г���������������ֵ%ע�⣡���ղ������ڷ�����
        rootCutVers = movedRootTooth2(rootCutIdx ,:);
        a_t1 = ismember(rootTooth.faces(:,1),rootCutIdx);
        aa_t1 = find(a_t1==1);
        a_t2 = ismember(rootTooth.faces(:,2),rootCutIdx);
        aa_t2 = find(a_t2==1);
        a_t3 = ismember(rootTooth.faces(:,3),rootCutIdx);
        aa_t3 = find(a_t3==1);
        [c_t1,~,~] = intersect(aa_t1,aa_t2);
        triIdxInCutRoot = intersect(c_t1,aa_t3);

        %�ұ߽�
        [raw_edges_list_t] = query_edges_list(rootTooth.faces(triIdxInCutRoot,:),'sorted');
        [~,~,iu_t] = unique(sort(raw_edges_list_t,2),'rows');
        [i3,i4] = histc(iu_t,unique(iu_t));
        lone_edges_idx_vect_t = i3(i4) == 1;
        rootBdryEdges = unique(raw_edges_list_t(lone_edges_idx_vect_t,:),'rows');
        
        
%       3.3 ȷ���и���������Ƭ���߽�ӳ�䵽�ϲ������е�����Ƭ���߽硣
        trisInCutRoot = rootTooth.faces(triIdxInCutRoot,:);
        rootCutTrisInMergeMesh = ones (size(trisInCutRoot));
        for j = 1:length(rootCutVers)
            nu = find( trisInCutRoot == rootCutIdx(j));
            rootCutTrisInMergeMesh(sub2ind(size(rootCutTrisInMergeMesh), nu)) = j+length(cutPatientCrownVers);
        end
        allTris = [f_b0; rootCutTrisInMergeMesh];
        allVers = [cutPatientCrownVers; rootCutVers];

       
        %�߽�㣺
        edge_p_t_local_list = unique([rootBdryEdges(:,1);rootBdryEdges(:,2)]);
        
        
        %�߽�����������е�λ��
        [C_t,edge_p_t_local,~] = intersect(rootCutIdx,edge_p_t_local_list);
  
        edge_p_t = rootCutVers(edge_p_t_local,:);
        edge_p_t_local = edge_p_t_local+length(cutPatientCrownVers);
        trisInCutRoot = zeros(size(rootBdryEdges));
        for i = 1:length(edge_p_t_local)
            nu = find(rootBdryEdges == edge_p_t_local_list(i));
            trisInCutRoot(nu) = edge_p_t_local(i);
        end
        
        
        
        
%% 4. �ϲ���������Ƭ

%5.20����
%1.���еı߽�㰴��һ����������
        [~,startEdgeIdx] = mindis(allVers(trisInCutRoot(:,1),:),allVers(Edge_b(1,1),:),1);
        
        %�ж�[edg_b(1,��),Edge_t(r,1)]de ������[Edge_t(r,��),edg_b(1,1)]�ķ����Ƿ���ͬ
        %Ϊ�˷�ֹ���������������Edge_t(r,1)����ĵڶ�����
        [r2,r3] = ind2sub(size(trisInCutRoot),find(trisInCutRoot == trisInCutRoot(startEdgeIdx,2)));
        rr = r2(ismember(r2,startEdgeIdx)==0);
        rm = r3(ismember(r2,startEdgeIdx)==0);
        tri11 = allVers(patientEdge(3,1),:) - allVers(patientEdge(1,1),:);
        tri12 = allVers(trisInCutRoot(startEdgeIdx,1),:) - allVers(patientEdge(3,1),:);
        tri21 = allVers(trisInCutRoot(rr,3-rm),:) - allVers(trisInCutRoot(startEdgeIdx,1),:);tri22 = allVers(patientEdge(1,1),:) - allVers(trisInCutRoot(rr,3-rm),:);
        s_b = sign(dot(allVers(patientEdge(1,1),:) -mean(edge_p_b),cross(tri11,tri12)));
        s_t = sign(dot(allVers(trisInCutRoot(startEdgeIdx,1),:) - mean(edge_p_t),cross(tri21,tri22)));
        if s_b ~= s_t
           edg_t = sotr_edge(trisInCutRoot,startEdgeIdx);
        else
            E_t = [trisInCutRoot(:,2),trisInCutRoot(:,1)];
            edg_t = sotr_edge(E_t,startEdgeIdx);
        end


        %�߽���ں�,��ϡ�裨�������ı���Ϊ��׼�����ܼ������ڣ������������������ĵ�
        
%   4.2 
        addTris = [];
        middle = [];
for j = 1:length(edg_t)
    pe1 = allVers(edg_t(j,1),:);
    p2 = allVers(edg_t(j,2),:);
    middle(j, :) = (pe1+p2)/2;
    [~,ro]=mindis(allVers(patientEdge(:,1),:),(pe1+p2)/2,1);
    roww(j) = ro;
    addTris = [addTris;[edg_t(j,:),patientEdge(ro,1)]];
end
 
OBJwriteVertices('middle.obj',middle);
        
        if roww(1)<length(patientEdge)/2
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
                addTris = [addTris;[patientEdge(roww(end):length(patientEdge),:),repmat(edg_t(1,1),length(patientEdge)-roww(end)+1,1)]];
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
        ff = [allTris;addTris];
        p = allVers;
        b = select_holes_and_boundary(allVers,ff);
        

%% 5  ���ο��Ʋ���
        index_crown_change = find(patientTooth.vertex(:,2)<(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+1 ...
                                & patientTooth.vertex(:,2)>=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-2);
        crownforchange = patientTooth.vertex(index_crown_change,:);
        
        index_tooth_change =  find(p(:,2)<(centerPatient(2) - ( patientYdir(1)*(p(:,1) ...
            - p_patient(1))+ patientYdir(3)*(p(:,3) - p_patient(3)))/ patientYdir(2))+1 ...
            &p(:,2)>=(p_patient(2) - ( patientYdir(1)*(p(:,1) ...
            - p_patient(1))+ patientYdir(3)*(p(:,3) - p_patient(3)))/ patientYdir(2))-2);
        pchange = p(index_tooth_change,:);
        
        
        %�������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�tooth_T����crownforchange��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
        p_T = p;
        bi_V = zeros(size(p_T));
        bi_bndtype = 'ext';
        BZ1 = zeros(size(p_T,1),3);
        reduction = 'no_flatten';masstype = 'voronoi';
        indices = 1:size(p_T,1);
        %�����ε�����
        exterior = indices(p(:,2)>(p_patient(2) - ( patientYdir(1)*(p(:,1) ...
            - p_patient(1))+ patientYdir(3)*(p(:,3) - p_patient(3)))/ patientYdir(2))+1.5 ...
            |p(:,2)<=(p_patient(2) - ( patientYdir(1)*(p(:,1) ...
            - p_patient(1))+ patientYdir(3)*(p(:,3) - p_patient(3)))/ patientYdir(2))-3);
    end
end
    


%% 6
[omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(p_T,1), ff, exterior);
[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(p,ff, bi_bndtype,masstype,reduction,omega,N0,N1);
%�ҵ��������ָ����ڱ��β�������ĵ㣬�������ϵĵ��������ڴ�
for i = 1:length(pchange)
   [minValue,startEdgeIdx]=mindis(crownforchange,pchange(i,:),1);
   minvalue(i) = minValue;  row(i) = startEdgeIdx;
   p_T(index_tooth_change(i),:) = crownforchange(startEdgeIdx,:);
end
bi_V = biharm_solve_with_factor( ...
    bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
    ff, p_T, omega, N0, N1, bi_bndtype, reduction,BZ1,p);
figure()
trisurf(ff,bi_V(:,1),bi_V(:,2),bi_V(:,3),'facecolor','c','edgecolor','b')
% axis image
 
writeOBJ('��������.obj',bi_V, ff)
 
