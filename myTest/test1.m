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
        rootTooth = dentalwithtooth1(i);         % �������ı�׼������
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
crownTooth = crown.model;
writeOBJ('��׼����.obj', crownTooth.vertex, crownTooth.face);
writeOBJ('��������׼��.obj', rootTooth.vertices, rootTooth.faces);        % �������������ı�׼����������Ƭ�Ƿ��ģ��Ƿ������⣿
crown_center = mean(crownTooth.vertex);
crown_ax_xline = getDirLine(crown_ax.x, crown_center);
crown_ax_yline = getDirLine(crown_ax.y, crown_center);
crown_ax_zline = getDirLine(crown_ax.z, crown_center);
OBJwriteVertices('��׼��x��.obj', crown_ax_xline);
OBJwriteVertices('��׼��y��.obj', crown_ax_yline);
OBJwriteVertices('��׼��z��.obj', crown_ax_zline);
crown_ax_dirs = [crown_ax.x; crown_ax.y; crown_ax.z];
OBJwriteVertices('��׼�����᷽������.obj', crown_ax_dirs);


if s ==1 || s ==2               % �����
    fdi = textread('FDIUpper__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
        disp( '�ò���û�д����� ');
        return;
       
    else
        namestr1 = ['toothUpper_',num2str(toothIdx-1),'.','obj'];
        patientTooth = Read_Obj(namestr1);                       % ���˵���������
        namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);                       % �����ߵ���
        axis = ReadObj('AXISUpper_.obj');
        axi = axis( (3*(toothIdx-1) + 1) : 3*toothIdx, :);       % ��ǰ���ݵ��������᷽������
        OBJwriteVertices('�������᷽������.obj', axi);
        OBJwriteVertices('����������.obj', gumline);
        
        
        % 1.2 ȷ���������ڵ��и��
        centerPatient = mean(patientTooth.vertex);     % ����������������
        ymax = max(gumline(:,2));             % �����ߵ�����y�������ֵ
        p_bingren =[centerPatient(1), ymax, centerPatient(3)];
        toothYdir = axi(2,:);                    % ��ǰ�������ݵ�y�᷽��
        index_bingren = find(patientTooth.vertex(:,2) <= (p_bingren(2) - (toothYdir(1)*(patientTooth.vertex(:,1) ...
                - p_bingren(1)) + toothYdir(3)*(patientTooth.vertex(:,3) - p_bingren(3)))/toothYdir(2))-2);    % ������Ϊ2
        cutPatientVers = patientTooth.vertex(index_bingren,:);      % ������������ϵ�£�ֻ����y�����ϸ���p_bingren��2mm�ĵ㡣
        
        disp(size(index_bingren));      % for debug
        
        
       % FOR DEBUG
        OBJwriteVertices('p_bingren.obj', p_bingren);
        moved_p_bingren = p_bingren + 2 * toothYdir;
        writeOBJ('������������.obj', patientTooth.vertex, patientTooth.face);
        OBJwriteVertices('1.2moved_p_bingren.obj', moved_p_bingren);
        OBJwriteVertices('1.2��ȥ�ײ��Ĳ������ڵ㼯.obj', cutPatientVers);
        bingren_xline = getDirLine(axi(1,:), centerPatient);
        bingren_yline = getDirLine(axi(2,:), centerPatient);
        bingren_zline = getDirLine(axi(3,:), centerPatient);
        OBJwriteVertices('������x��.obj', bingren_xline);
        OBJwriteVertices('������y��.obj', bingren_yline);
        OBJwriteVertices('������z��.obj', bingren_zline);
        
        centerCrown = mean(crownTooth.vertex);     % ��׼��������
        crownYdir = crown_ax.y;
        index_crown = find(crownTooth.vertex(:,2)<=(centerCrown(2) - (crownYdir(1)*(crownTooth.vertex(:,1) ...
                    - centerCrown(1))+crownYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/crownYdir(2)));
        v_crown_biaozhun = abs(centerCrown(2) - (crownYdir(1)*(crownTooth.vertex(index_crown,1)-centerCrown(1))...
                     +crownYdir(3)*(crownTooth.vertex(index_crown,3)-centerCrown(3)))/crownYdir(2));
        in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
        point_crown =  crownTooth.vertex(index_crown(in_crown),:);
        
       %  for debug
        OBJwriteVertices('point_crown.obj', point_crown);
        
       % 1.3 ȷ��������׼���и��
        centerRoot = mean(rootTooth.vertices);
        index_root = find(rootTooth.vertices(:,2)<=(centerRoot(2) - (crownYdir(1)*(rootTooth.vertices(:,1) ...
                         - centerRoot(1))+crownYdir(3)*(rootTooth.vertices(:,3) ...
                         - centerRoot(3)))/crownYdir(2)));       % ������������ϵ�£�y��߶ȸ������ĵĵ㡣
        vz_tooth = abs(centerRoot(2) - (crownYdir(1)*(rootTooth.vertices(index_root,1)-centerRoot(1))...
                     +crownYdir(3)*(rootTooth.vertices(index_root,3)-centerRoot(3)))/crownYdir(2));
        in_tooth = find(vz_tooth == max(vz_tooth));
        point_tooth = rootTooth.vertices(index_root(in_tooth),:);
        centcrownintooth = centerCrown + point_tooth - point_crown;  % �к��׼���ı�Ե����
        
        % for debug
        OBJwriteVertices('centerRoot.obj', centerRoot);  
        OBJwriteVertices('point_tooth.obj', point_tooth);  
        OBJwriteVertices('�к��׼���ı�Ե����.obj', centcrownintooth);  
        
        
       % 2. ���롪�����ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж��룬������ղ������ڵ��и�λ�ý����и�
        R = inv([crown_ax.x; crown_ax.y; crown_ax.z]) * axi;   % ��ת
        C = centcrownintooth * R;         
        tooth_T = (rootTooth.vertices) *R + repmat((centerPatient - C),length(rootTooth.vertices),1);
        
        %       2.1 ��תƽ�ƴ�����׼��
        index_1 = find(patientTooth.vertex(:,2)>=(centerPatient(2) - (toothYdir(1)*(patientTooth.vertex(:,1) - p_bingren(1))+toothYdir(3)*(patientTooth.vertex(:,3) - p_bingren(3)))/toothYdir(2))...
                 & patientTooth.vertex(:,2)<=(p_bingren(2) - (toothYdir(1)*(patientTooth.vertex(:,1) - ...
                 p_bingren(1))+toothYdir(3)*(patientTooth.vertex(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        patientMergeVers = patientTooth.vertex(index_1,:);
        
        index_2 = find(tooth_T(:,2)>=(centerPatient(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2))...
                & tooth_T(:,2)<=(p_bingren(2) - (toothYdir(1)*(tooth_T(:,1) - p_bingren(1))+toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/toothYdir(2))+0.5);
        rootMergeVers = tooth_T(index_2,:);
        
        % for debug
         OBJwriteVertices('2.1��תƽ�ƺ�Ĵ�����׼��.obj', tooth_T);
         OBJwriteVertices('2.1�����������ں�����ĵ�.obj', patientMergeVers);
         OBJwriteVertices('2.1������׼���ں�����ĵ�.obj', rootMergeVers);
        
        %       2.2 ȷ��������׼�����и�֡�
        for i = 1:length(rootMergeVers)
           [minValue,r]=mindis(patientMergeVers, rootMergeVers(i,:),1);
           minvalue(i) = minValue;  row(i) = r;
        end
        point11 = patientMergeVers(row,:);
        [R,t,BRt,e,~,~] = icp(point11,rootMergeVers);
        tooth_T = bsxfun(@plus,tooth_T*R,t);        
        index_tooth =  find(tooth_T(:,2)>(p_bingren(2) - ( toothYdir(1)*(tooth_T(:,1) ...
                         - p_bingren(1))+ toothYdir(3)*(tooth_T(:,3) - p_bingren(3)))/ toothYdir(2))-2);  %�г���������������ֵ%ע�⣡���ղ������ڷ�����
        rootCutVers = tooth_T(index_tooth ,:);
      

        %       2.3 ���ڱ��ο��Ʋ���
        index_crown_change = find(patientTooth.vertex(:,2)>(centerPatient(2) - (toothYdir(1)*(patientTooth.vertex(:,1) - p_bingren(1))+toothYdir(3)*(patientTooth.vertex(:,3) - p_bingren(3)))/toothYdir(2))-1 ...
                                & patientTooth.vertex(:,2)<=(p_bingren(2) - (toothYdir(1)*(patientTooth.vertex(:,1) - p_bingren(1))+toothYdir(3)*(patientTooth.vertex(:,3) - p_bingren(3)))/toothYdir(2))+2);
        crownforchange = patientTooth.vertex(index_crown_change,:);
        
       % for debug
        OBJwriteVertices('�г�������������.obj', rootCutVers);

        %       2.5 ���ں������ϲ���һ���������񣿣�����
        mergedToothVers = [cutPatientVers; rootCutVers];
        [mergedToothFace]=MyCrustOpen(mergedToothVers);
        % for debug
        writeOBJ('�ϲ�����.Obj', mergedToothVers, mergedToothFace);             % ������ ����Ƭ�������

       
        %       2.6 ���񲹶�
        mergedToothBdr = select_holes_and_boundary(mergedToothVers, mergedToothFace);
        ff = fill_mesh_holes(mergedToothVers, mergedToothFace, mergedToothBdr,'closed',200);
        ff = double(ff);
        index_tooth_change =  find(mergedToothVers(:,2)>(centerPatient(2) - ( toothYdir(1)*(mergedToothVers(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedToothVers(:,3) - p_bingren(3)))/ toothYdir(2))-1 ...
            &mergedToothVers(:,2)<=(p_bingren(2) - ( toothYdir(1)*(mergedToothVers(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedToothVers(:,3) - p_bingren(3)))/ toothYdir(2))+1);
        pchange = mergedToothVers(index_tooth_change,:);
        
        
        %       2.7 �������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�tooth_T����crownforchange��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
        mergedToothVers_copy = mergedToothVers;
        BZ1 = zeros(size(mergedToothVers_copy,1),3);
        indices = 1:size(mergedToothVers_copy,1);
        
        %       2.8 �����ε�����
        exterior = indices(mergedToothVers(:,2)<(p_bingren(2) - ( toothYdir(1)*(mergedToothVers(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedToothVers(:,3) - p_bingren(3)))/ toothYdir(2))-0.6...
            |mergedToothVers(:,2)>=(p_bingren(2) - ( toothYdir(1)*(mergedToothVers(:,1) ...
            - p_bingren(1))+ toothYdir(3)*(mergedToothVers(:,3) - p_bingren(3)))/ toothYdir(2))+3);
    end
else                             % �����
    % ��������������
end



[Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(mergedToothVers_copy,1), ff, exterior);


[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(mergedToothVers,ff, 'ext', 'voronoi', 'no_flatten',Omega,N0,N1);


%�ҵ��������ָ����ڱ��β�������ĵ㣬�������ϵĵ��������ڴ�
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

 
