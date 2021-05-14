clc
clear all
functionname='testforfusion.m'; 
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
% load('up_bingren.mat');
% load('low_bingren.mat');

%��ͬ��FDI���ó���
prompt  = '������FDI��� \n';
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
s = fix(x/10);  %ȡ��
g = mod(x,10);%ȡ��
if s ==1 || s ==2
%     addpath([functiondir 'data2\53120\u'])
%     fdi = textread('FDI__.txt');
    fdi = textread('FDIUpper__.dxt');
    n = find (fdi == x);
    if isempty(n)
        disp( '�ò���û�д����� ');
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
        %����
        c = mean(bingren.vertex);y = max(gumline(:,2));
        p_bingren =[c(1),y(1),c(3)];n_bingren = axi(2,:);v_bingren = bingren.vertex;
        index_bingren = find(v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-2);
        point_bingren = v_bingren(index_bingren,:);
%         figure()
%         plot3(point_bingren(:,1),point_bingren(:,2),point_bingren(:,3),'g*')
%         hold on

        %��׼����
        face_crown_biaozhun = crown.model.face;vertex_crown_biaozhun = crown.model.vertex;
        p_crown_biaozhun = mean(vertex_crown_biaozhun);
        n_crown_biaozhun = crown_ax.y;
        index_crown_biaozhun = find(vertex_crown_biaozhun(:,2)<=(p_crown_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_crown_biaozhun(:,1) ...
                    - p_crown_biaozhun(1))+n_crown_biaozhun(3)*(vertex_crown_biaozhun(:,3) - p_crown_biaozhun(3)))/n_crown_biaozhun(2)));
        v_crown_biaozhun = abs(p_crown_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_crown_biaozhun(index_crown_biaozhun,1)-p_crown_biaozhun(1))...
                     +n_crown_biaozhun(3)*(vertex_crown_biaozhun(index_crown_biaozhun,3)-p_crown_biaozhun(3)))/n_crown_biaozhun(2));
        in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
        point_crown =  vertex_crown_biaozhun(index_crown_biaozhun(in_crown),:);

        %��׼����
        vertex_tooth_biaozhun = root.vertices;
        f = root.faces;
        p_tooth_biaozhun = mean(vertex_tooth_biaozhun);
        index_tooth_biaozhun = find(vertex_tooth_biaozhun(:,2)<=(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(:,1) ...
                         - p_tooth_biaozhun(1))+n_crown_biaozhun(3)*(vertex_tooth_biaozhun(:,3) - p_tooth_biaozhun(3)))/n_crown_biaozhun(2)));
        vz_tooth = abs(p_tooth_biaozhun(2) - (n_crown_biaozhun(1)*(vertex_tooth_biaozhun(index_tooth_biaozhun,1)-p_tooth_biaozhun(1))...
                   +n_crown_biaozhun(3)*(vertex_tooth_biaozhun(index_tooth_biaozhun,3)-p_tooth_biaozhun(3)))/n_crown_biaozhun(2));
        in_tooth = find(vz_tooth == max(vz_tooth));
        point_tooth = vertex_tooth_biaozhun(index_tooth_biaozhun(in_tooth),:);
        centcrownintooth = p_crown_biaozhun+point_tooth-point_crown;%�к��׼���ı�Ե����
        %%���ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж��룬������ղ������ڵ��и�λ�ý����и�
        R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*axi;
        C=centcrownintooth*R;
        tooth_T = (vertex_tooth_biaozhun) *R + repmat((c - C),length(vertex_tooth_biaozhun),1);
        
        %����仯���ݵ���̬
        index_1 = find(v_bingren(:,2)>=(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))...
                & v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+0.5);
        point1 = v_bingren(index_1,:);
        index_2 = find(tooth_T(:,2)>=(c(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2))...
                & tooth_T(:,2)<=(p_bingren(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2))+0.5);
        point2 = tooth_T(index_2,:);
        %ȷ�����β���
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
                         - p_bingren(1))+ n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/ n_bingren(2))-2);%�г���������������ֵ%ע�⣡���ղ������ڷ�����
        tooth = tooth_T(index_tooth ,:);
  
%         plot3(tooth(:,1),tooth(:,2),tooth(:,3),'r*')

%         tooth = vertex_tooth_biaozhun(index_tooth,:);
%         R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*bingren.orientation;
%         C=centcrownintooth*R;
%         tooth_T = (tooth) *R + repmat((p_bingren - C),length(tooth),1);

        %���ڱ��ο��Ʋ���
        index_crown_change = find(v_bingren(:,2)>(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-1 ...
                                & v_bingren(:,2)<=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+2);
        crownforchange = v_bingren(index_crown_change,:);
        %���������ֱ��γ�����
        [tooth_t]=MyCrustOpen(point_bingren);
        figure()
        trimesh(tooth_t,point_bingren(:,1),point_bingren(:,2),point_bingren(:,3))
        hold on
        [tooth_t2]=MyCrustOpen(tooth); 
        trimesh(tooth_t2,tooth(:,1),tooth(:,2),tooth(:,3))
        hold off
        %���ں������γ�����
        p = [point_bingren;tooth];
        [t]=MyCrustOpen(p);
        figure()
        trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
        %���񲹶�
        b = select_holes_and_boundary(p,t);
        ff = fill_mesh_holes(p,t,b,'closed',200);
        ff = double(ff);
        index_tooth_change =  find(p(:,2)>(c(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-1 ...
            &p(:,2)<=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+1);
        pchange = p(index_tooth_change,:);
        %�������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�tooth_T����crownforchange��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
        p_T = p;
        bi_V = zeros(size(p_T));
        bi_bndtype = 'ext';
        BZ1 = zeros(size(p_T,1),3);
        reduction = 'no_flatten';masstype = 'voronoi';
        indices = 1:size(p_T,1);
        %�����ε�����
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
       disp( '�ò���û�д����� ');
       return;
    else
        namestr1 = ['toothLower_',num2str(n-1),'.','obj'];
        bingren = Read_Obj(namestr1);
        namestr2 = ['gumlineLower_',num2str(n-1),'.','obj'];
        gumline =  ReadObj(namestr2);
        axis = ReadObj('AXISLower_.obj');
        axi = axis(3*(n-1)+1:3*n,:);
    
         %����
        c = mean(bingren.vertex);y = min(gumline(:,2));
        p_bingren =[c(1),y(1),c(3)];n_bingren = axi(2,:);v_bingren = bingren.vertex;
        index_bingren = find(v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+1);
        point_bingren = v_bingren(index_bingren,:);
%         figure()
%         plot3(point_bingren(:,1),point_bingren(:,2),point_bingren(:,3),'g*')
%         hold on

        %��׼����
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

        %��׼����
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
%         centcrownintooth = p_crown_biaozhun+point_tooth-point_crown;%�к��׼���ı�Ե����
        vertex_tooth_biaozhun = vertex_tooth_biaozhun-point_tooth+point_crown;
        %%���ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж��룬������ղ������ڵ��и�λ�ý����и�
        R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*axi;
        C=p_crown_biaozhun*R;
        tooth_T = (vertex_tooth_biaozhun) *R + repmat((c - C),length(vertex_tooth_biaozhun),1);
        figure()
        trisurf(f,tooth_T(:,1),tooth_T(:,2),tooth_T(:,3),'facecolor','c','edgecolor','b')
        hold on
        trisurf(bingren.face,v_bingren(:,1),v_bingren(:,2),v_bingren(:,3),'facecolor','c','edgecolor','r')
        
        %����仯���ݵ���̬
        index_1 = find(v_bingren(:,2)<=(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2)) ...
                & v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+0.5);
        point1 = v_bingren(index_1,:);
        index_2 = find(tooth_T(:,2)<=(c(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2)) ...
                & tooth_T(:,2)>=(p_bingren(2) - (n_bingren(1)*(tooth_T(:,1) - p_bingren(1))+n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/n_bingren(2))+0.5);
        point2 = tooth_T(index_2,:);
        %ȷ�����β���
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
                         - p_bingren(1))+ n_bingren(3)*(tooth_T(:,3) - p_bingren(3)))/ n_bingren(2))+1);%�г���������������ֵ%ע�⣡���ղ������ڷ�����
        tooth = tooth_T(index_tooth ,:);
        
        
        
        
        
        
%         plot3(tooth(:,1),tooth(:,2),tooth(:,3),'r*')

%         tooth = vertex_tooth_biaozhun(index_tooth,:);
%         R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*bingren.orientation;
%         C=centcrownintooth*R;
%         tooth_T = (tooth) *R + repmat((p_bingren - C),length(tooth),1);

        %���ڱ��ο��Ʋ���
        index_crown_change = find(v_bingren(:,2)<(c(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))+1 ...
                                & v_bingren(:,2)>=(p_bingren(2) - (n_bingren(1)*(v_bingren(:,1) - p_bingren(1))+n_bingren(3)*(v_bingren(:,3) - p_bingren(3)))/n_bingren(2))-1);
        crownforchange = v_bingren(index_crown_change,:);
        %���������ֱ��γ�����
        [tooth_t]=MyCrustOpen(point_bingren);
        figure()
        trimesh(tooth_t,point_bingren(:,1),point_bingren(:,2),point_bingren(:,3))
        hold on
        [tooth_t2]=MyCrustOpen(tooth); 
        trimesh(tooth_t2,tooth(:,1),tooth(:,2),tooth(:,3))
        hold off
        %���ں������γ�����
        p = [point_bingren;tooth];
        [t]=MyCrustOpen(p);
        figure()
        trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')
        %���񲹶�
        b = select_holes_and_boundary(p,t);
        ff = fill_mesh_holes(p,t,b,'closed',200);
        ff = double(ff);
        index_tooth_change =  find(p(:,2)<(c(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+1 ...
            &p(:,2)>=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-1);
        pchange = p(index_tooth_change,:);
        %�������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�tooth_T����crownforchange��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
        p_T = p;
        bi_V = zeros(size(p_T));
        bi_bndtype = 'ext';
        BZ1 = zeros(size(p_T,1),3);
        reduction = 'no_flatten';masstype = 'voronoi';
        indices = 1:size(p_T,1);
        %�����ε�����
        exterior = indices(p(:,2)>(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))+1 ...
            |p(:,2)<=(p_bingren(2) - ( n_bingren(1)*(p(:,1) ...
            - p_bingren(1))+ n_bingren(3)*(p(:,3) - p_bingren(3)))/ n_bingren(2))-3);
    end
    
    
    
end


[Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(p_T,1), ff, exterior);
[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(p,ff, bi_bndtype,masstype,reduction,Omega,N0,N1);
%�ҵ��������ָ����ڱ��β�������ĵ㣬�������ϵĵ��������ڴ�
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
%     %�г�����
%     
%     
%     %�г���������,���ų��᷽�����볤�ᴹֱ��ƽ���ϣ��е�ƽ���ϵĵ�ȷ�����⣬��ʱ�����ڵ����ĵ㣬��������
%     
%     quiver3(up(i+1).center(1),up(i+1).center(2),up(i+1).center(3),up(i+1).orientation(3,1),up(i+1).orientation(3,2),up(i+1).orientation(3,3),10); 
%     
%     %������루����һ���̶�����ᣬ����ƽ������ĵ㣩
%     
% end

