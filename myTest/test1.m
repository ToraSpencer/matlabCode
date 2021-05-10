clc
close all;
clear all;
% functionname='test1.m'; 
% functiondir=which(functionname);
% functiondir=functiondir(1:end-length(functionname));
% addpath([functiondir 'mesh process'])
addpath(['mesh process'])

disp('start')
%% 1. ���ر�׼��
load('dental_crown.mat');                       % ��׼����ֻ������û������
load('dentalmodelwithroot0.1forlow.mat');       % ��׼����������
load('axisofdental_crowm');                     %            
 

%       1.1 ��ͬ��FDI���ó���
x = 11;             % ȡ11����������
for i = 1:28
    if  dentalwithtooth1(i).ID == x
        root = dentalwithtooth1(i);         % �������ı�׼������
    end
end
 
for j = 1:28
    if (upax(j).fid == x)
        crown_ax = upax(j);
    end
end

for k = 1:28
    if (dental_crown(k).fid == x)
        crown = dental_crown(k);            % ֻ�����ڵı�׼�����񣿣���
    end
end


s = fix(x/10);  %ȡ��
g = mod(x,10); %ȡ��


%% for debug
writeOBJ('��׼����.obj', crown.model.vertex, crown.model.face);
writeOBJ('��������׼��.obj', root.vertices, root.faces);        % �������������ı�׼����������Ƭ�Ƿ��ģ��Ƿ������⣿
crown_center = mean(crown.model.vertex);
crown_ax_xline = getDirLine(crown_ax.x, crown_center);
crown_ax_yline = getDirLine(crown_ax.y, crown_center);
crown_ax_zline = getDirLine(crown_ax.z, crown_center);
OBJwriteVertices('��׼��x��.obj', crown_ax_xline);
OBJwriteVertices('��׼��y��.obj', crown_ax_yline);
OBJwriteVertices('��׼��z��.obj', crown_ax_zline);
crown_ax_dirs = [crown_ax.x; crown_ax.y; crown_ax.z];
OBJwriteVertices('��׼�����᷽������.obj', crown_ax_dirs);




if s ==1 || s ==2
    fdi = textread('FDIUpper__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
        disp( '�ò���û�д����� ');
        return;
       
    else
        namestr1 = ['toothUpper_',num2str(toothIdx-1),'.','obj'];
        bingren = Read_Obj(namestr1);                       % ���˵���������
        namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);                       % �����ߵ���
        axis = ReadObj('AXISUpper_.obj');
        axi = axis( (3*(toothIdx-1) + 1) : 3*toothIdx, :);       % ��ǰ���ݵ��������᷽������
        OBJwriteVertices('�������᷽������.obj', axi);
        OBJwriteVertices('����������.obj', gumline);
        %       1.2 ����
        bingren_center = mean(bingren.vertex);     % ����������������
        ymax = max(gumline(:,2));             % �����ߵ�����y�������ֵ
        p_bingren =[bingren_center(1), ymax, bingren_center(3)];
        toothYdir = axi(2,:);                    % ��ǰ�������ݵ�y�᷽��
        bingren_vers = bingren.vertex;
        index_bingren = find(bingren_vers(:,2) <= (p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) ...
                - p_bingren(1)) + toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))-2);    % ������Ϊ2
        point_bingren = bingren_vers(index_bingren,:);      % ������������ϵ�£�ֻ����y�����ϸ���p_bingren��2mm�ĵ㡣
        
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
        
        writeOBJ('��������.obj', bingren.vertex, bingren.face);
        OBJwriteVertices('moved_p_bingren.obj', moved_p_bingren);
        OBJwriteVertices('��ȥ�ײ��Ĳ������ڵ㼯.obj', point_bingren);
        
        bingren_xline = getDirLine(axi(1,:), bingren_center);
        bingren_yline = getDirLine(axi(2,:), bingren_center);
        bingren_zline = getDirLine(axi(3,:), bingren_center);
        OBJwriteVertices('������x��.obj', bingren_xline);
        OBJwriteVertices('������y��.obj', bingren_yline);
        OBJwriteVertices('������z��.obj', bingren_zline);
        
       %%
        
        %       ��׼����
        face_crown_biaozhun = crown.model.face;
        vertex_crown_biaozhun = crown.model.vertex;
        crown_biaozhun_center = mean(vertex_crown_biaozhun);     % ��׼��������
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
       
        % ��׼����
        rootVers = root.vertices;
        root_center = mean(rootVers);
        index_root = find(rootVers(:,2)<=(root_center(2) - (biaozhunYdir(1)*(rootVers(:,1) ...
                         - root_center(1))+biaozhunYdir(3)*(rootVers(:,3) ...
                         - root_center(3)))/biaozhunYdir(2)));       % ������������ϵ�£�y��߶ȸ������ĵĵ㡣
        vz_tooth = abs(root_center(2) - (biaozhunYdir(1)*(rootVers(index_root,1)-root_center(1))...
                     +biaozhunYdir(3)*(rootVers(index_root,3)-root_center(3)))/biaozhunYdir(2));
        in_tooth = find(vz_tooth == max(vz_tooth));
        point_tooth = rootVers(index_root(in_tooth),:);
        centcrownintooth = crown_biaozhun_center + point_tooth - point_crown;  % �к��׼���ı�Ե����
        
        %% for debug
        chosenRootVers = rootVers(index_root, :);
        OBJwriteVertices('root_center.obj', root_center);  
        OBJwriteVertices('chosenRootVers.obj', chosenRootVers);      % ������������ϵ�£�y��߶ȸ������ĵĵ㡣
        OBJwriteVertices('point_tooth.obj', point_tooth);  
        OBJwriteVertices('�к��׼���ı�Ե����.obj', centcrownintooth);  
        %%
        
       %% 2. ���롪�����ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж��룬������ղ������ڵ��и�λ�ý����и�
        R = inv([crown_ax.x; crown_ax.y; crown_ax.z]) * axi;   % ��ת
        C = centcrownintooth * R;         
        tooth_T = (rootVers) *R + repmat((bingren_center - C),length(rootVers),1);
        
        %       2.1 ����仯���ݵ���̬
        index_1 = find(bingren_vers(:,2)>=(bingren_center(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))...
                 & bingren_vers(:,2)<=(p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) - ...
                 p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        point1 = bingren_vers(index_1,:);
        
        index_2 = find(tooth_T(:,2)>=(bingren_center(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2))...
                & tooth_T(:,2)<=(p_bingren(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        point2 = tooth_T(index_2,:);
        
        %% for debug
         OBJwriteVertices('��תƽ�ƺ�Ĵ�����׼��.obj', tooth_T);
         OBJwriteVertices('point1.obj', point1);
         OBJwriteVertices('point2.obj', point2);
        
       %%
        %       2.2 ȷ�����β���
        for i = 1:length(point2)
           [minValue,r]=mindis(point1,point2(i,:),1);
           minvalue(i) = minValue;  row(i) = r;
        end
        point11 = point1(row,:);
        [R,t,BRt,e,~,~] = icp(point11,point2);
        tooth_T = bsxfun(@plus,tooth_T*R,t);        
        index_tooth =  find(tooth_T(:,2)>(p_bingren(2) - ( toothYdir(1)*(tooth_T(:,1) ...
                         - p_bingren(1))+ toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/ toothYdir(2))-2);%�г���������������ֵ%ע�⣡���ղ������ڷ�����
        rootCutVers = tooth_T(index_tooth ,:);
        

        

        %       2.3 ���ڱ��ο��Ʋ���
        index_crown_change = find(bingren_vers(:,2)>(bingren_center(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))-1 ...
                                & bingren_vers(:,2)<=(p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))+2);
        crownforchange = bingren_vers(index_crown_change,:);
        
       %% for debug
        OBJwriteVertices('�г�������������.obj', rootCutVers);
        figure()
        trimesh(MyCrustOpen(point_bingren),point_bingren(:,1),point_bingren(:,2),point_bingren(:,3))
        hold on
        trimesh(MyCrustOpen(rootCutVers),rootCutVers(:,1),rootCutVers(:,2),rootCutVers(:,3))
        hold off
       %%
        
        %       2.5 ���ں������ϲ���һ���������񣿣�����
        mergedTooth = [point_bingren; rootCutVers];
        [mergedTooth_face]=MyCrustOpen(mergedTooth);
        figure()
        trisurf(mergedTooth_face,mergedTooth(:,1),mergedTooth(:,2),mergedTooth(:,3),'facecolor','c','edgecolor','b');
        
        %% for debug
        writeOBJ('�ϲ�����.Obj', mergedTooth, mergedTooth_face);             % ������ ����Ƭ�������
        %%
       
        %       2.6 ���񲹶�
        b = select_holes_and_boundary(mergedTooth, mergedTooth_face);
        ff = fill_mesh_holes(mergedTooth, mergedTooth_face, b,'closed',200);
        ff = double(ff);
        index_tooth_change =  find(mergedTooth(:,2)>(bingren_center(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))-1 ...
            &mergedTooth(:,2)<=(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))+1);
        pchange = mergedTooth(index_tooth_change,:);
        
        
        %       2.7 �������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�tooth_T����crownforchange��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
        p_T = mergedTooth;
        bi_V = zeros(size(p_T));
        bi_bndtype = 'ext';
        BZ1 = zeros(size(p_T,1),3);
        reduction = 'no_flatten';masstype = 'voronoi';
        indices = 1:size(p_T,1);
        
        %       2.8 �����ε�����
        exterior = indices(mergedTooth(:,2)<(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))-0.6...
            |mergedTooth(:,2)>=(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))+3);
    end
else
    fdi = textread('FDILower__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
       disp( '�ò���û�д����� ');
       return;
    else
        namestr1 = ['toothLower_',num2str(toothIdx-1),'.','obj'];
        bingren = Read_Obj(namestr1);
        namestr2 = ['gumlineLower_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);
        axis = ReadObj('AXISLower_.obj');
        axi = axis(3*(toothIdx-1)+1:3*toothIdx,:);
    
         %����
        bingren_center = mean(bingren.vertex);
        y = min(gumline(:,2));
        p_bingren =[bingren_center(1),y(1),bingren_center(3)];toothYdir = axi(2,:);bingren_vers = bingren.vertex;
        index_bingren = find(bingren_vers(:,2)>=(p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))+1);
        point_bingren = bingren_vers(index_bingren,:);
        
        %��׼����
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
 

        %��׼����
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
        
        
        %%���ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж��룬������ղ������ڵ��и�λ�ý����и�
        R = inv([crown_ax.x;crown_ax.y;crown_ax.z])*axi;
        C=crown_biaozhun_center*R;
        tooth_T = (rootVers) *R + repmat((bingren_center - C),length(rootVers),1);
        figure()
        trisurf(f,tooth_T(:,1),tooth_T(:,2),tooth_T(:,3),'facecolor','c','edgecolor','b')
        hold on
        trisurf(bingren.face,bingren_vers(:,1),bingren_vers(:,2),bingren_vers(:,3),'facecolor','c','edgecolor','r')
        
        %����仯���ݵ���̬
        index_1 = find(bingren_vers(:,2)<=(bingren_center(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2)) ...
                & bingren_vers(:,2)>=(p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        point1 = bingren_vers(index_1,:);
        index_2 = find(tooth_T(:,2)<=(bingren_center(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2)) ...
                & tooth_T(:,2)>=(p_bingren(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        point2 = tooth_T(index_2,:);
        
        
        %ȷ�����β���
        for i = 1:length(point2)
           [minValue,r]=mindis(point1,point2(i,:),1);
           minvalue(i) = minValue;  row(i) = r;
        end
        point11 = point1(row,:);
        [R,t,BRt,e,~,~] = icp(point11,point2);
        tooth_T = bsxfun(@plus,tooth_T,t);
        
        index_tooth =  find(tooth_T(:,2)<(p_bingren(2) - ( toothYdir(1)*(tooth_T(:,1) ...
                         - p_bingren(1))+ toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/ toothYdir(2))+1);%�г���������������ֵ%ע�⣡���ղ������ڷ�����
        rootCutVers = tooth_T(index_tooth ,:);
        
 
        %���ڱ��ο��Ʋ���
        index_crown_change = find(bingren_vers(:,2)<(bingren_center(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))+1 ...
                                & bingren_vers(:,2)>=(p_bingren(2) - (toothYdir(1)*(bingren_vers(:,1) - p_bingren(1))+toothYdir(3)*(bingren_vers(:,3) - p_bingren(3)))/toothYdir(2))-1);
        crownforchange = bingren_vers(index_crown_change,:);
        
        %���������ֱ��γ�����
        [tooth_t]=MyCrustOpen(point_bingren);
        figure()
        trimesh(tooth_t,point_bingren(:,1),point_bingren(:,2),point_bingren(:,3))
        hold on
        [tooth_t2]=MyCrustOpen(rootCutVers); 
        trimesh(tooth_t2,rootCutVers(:,1),rootCutVers(:,2),rootCutVers(:,3))
        hold off
       
        %���ں������γ�����
        mergedTooth = [point_bingren;rootCutVers];
        [t]=MyCrustOpen(mergedTooth);
        figure()
        trisurf(t,mergedTooth(:,1),mergedTooth(:,2),mergedTooth(:,3),'facecolor','c','edgecolor','b')
        
        %���񲹶�
        b = select_holes_and_boundary(mergedTooth,t);   % ��������ǰ���е�����ᱨ��������matlab2015�����ݵ����⡣
        ff = fill_mesh_holes(mergedTooth,t,b,'closed',200);
        ff = double(ff);
        index_tooth_change =  find(mergedTooth(:,2)<(bingren_center(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))+1 ...
            &mergedTooth(:,2)>=(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))-1);
        pchange = mergedTooth(index_tooth_change,:);
        
        %�������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�tooth_T����crownforchange��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
        p_T = mergedTooth;
        bi_V = zeros(size(p_T));
        bi_bndtype = 'ext';
        BZ1 = zeros(size(p_T,1),3);
        reduction = 'no_flatten';masstype = 'voronoi';
        indices = 1:size(p_T,1);
        
        %�����ε�����
        exterior = indices(mergedTooth(:,2)>(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))+1 ...
            |mergedTooth(:,2)<=(p_bingren(2) - ( toothYdir(1)*(mergedTooth(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedTooth(:,3) - p_bingren(3)))/ toothYdir(2))-3);
    end
    
    
    
end


[Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(p_T,1), ff, exterior);
[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(mergedTooth,ff, bi_bndtype,masstype,reduction,Omega,N0,N1);
%�ҵ��������ָ����ڱ��β�������ĵ㣬�������ϵĵ��������ڴ�
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

 
