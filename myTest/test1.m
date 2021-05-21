clc
close all;
clear all;
% functionname='test1.m'; 
% functiondir=which(functionname);
% functiondir=functiondir(1:end-length(functionname));
% addpath([functiondir 'mesh process'])
addpath('mesh process');
addpath('mixFE');

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
    fdi = textread('FDIUpper__.dxt');
    toothIdx = find (fdi == x);
    
    
namestr1 = ['toothUpper_',num2str(toothIdx-1),'.','obj'];
patientTooth = Read_Obj(namestr1);                       % ���˵���������
namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];
gumline =  ReadObj(namestr2);                       % �����ߵ���
axis = ReadObj('AXISUpper_.obj');
axi = axis( (3*(toothIdx-1) + 1) : 3*toothIdx, :);       % ��ǰ���ݵ��������᷽������
OBJwriteVertices('�������᷽������.obj', axi);
OBJwriteVertices('����������.obj', gumline);


%       1.2 ȷ���������ڵ��и��
centerPatient = mean(patientTooth.vertex);     % ����������������
ymax = max(gumline(:,2));             % �����ߵ�����y�������ֵ
p_patient =[centerPatient(1), ymax, centerPatient(3)];
patientYdir = axi(2,:);                    % ��ǰ�������ݵ�y�᷽��
index_bingren = find(patientTooth.vertex(:,2) <= (p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) ...
        - p_patient(1)) + patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-2);    % ������Ϊ2
cutPatientVers = patientTooth.vertex(index_bingren,:);      % ������������ϵ�£�ֻ����y�����ϸ���p_bingren��2mm�ĵ㡣

disp(size(index_bingren));      % for debug


% FOR DEBUG
OBJwriteVertices('p_bingren.obj', p_patient);
moved_p_bingren = p_patient + 2 * patientYdir;
writeOBJ('������������.obj', patientTooth.vertex, patientTooth.face);
OBJwriteVertices('moved_p_bingren.obj', moved_p_bingren);
OBJwriteVertices('�ϲ�����Ĳ������ڲ��ֵ㼯.obj', cutPatientVers);
bingren_xline = getDirLine(axi(1,:), centerPatient);
bingren_yline = getDirLine(axi(2,:), centerPatient);
bingren_zline = getDirLine(axi(3,:), centerPatient);
OBJwriteVertices('������x��.obj', bingren_xline);
OBJwriteVertices('������y��.obj', bingren_yline);
OBJwriteVertices('������z��.obj', bingren_zline);
OBJwriteVertices('����������.obj', centerPatient);

centerCrown = mean(crownTooth.vertex);     % ��׼��������
rootYdir = crown_ax.y;
index_crown = find(crownTooth.vertex(:,2)<=(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(:,1) ...
            - centerCrown(1))+rootYdir(3)*(crownTooth.vertex(:,3) - centerCrown(3)))/rootYdir(2)));
v_crown_biaozhun = abs(centerCrown(2) - (rootYdir(1)*(crownTooth.vertex(index_crown,1)-centerCrown(1))...
             +rootYdir(3)*(crownTooth.vertex(index_crown,3)-centerCrown(3)))/rootYdir(2));
in_crown = find(v_crown_biaozhun == max(v_crown_biaozhun));
point_crown =  crownTooth.vertex(index_crown(in_crown),:);

%  for debug
OBJwriteVertices('point_crown.obj', point_crown);

%        1.3 ȷ��������׼���и��
centerRoot = mean(rootTooth.vertices);
index_root_higher = find(rootTooth.vertices(:,2)<=(centerRoot(2) - (rootYdir(1)*(rootTooth.vertices(:,1) ...
                 - centerRoot(1))+rootYdir(3)*(rootTooth.vertices(:,3) ...
                 - centerRoot(3)))/rootYdir(2)));       % ������������ϵ�£�y��߶ȸ������ĵĵ㡣
rootZvalue = abs(centerRoot(2) - (rootYdir(1)*(rootTooth.vertices(index_root_higher,1)-centerRoot(1))...
             +rootYdir(3)*(rootTooth.vertices(index_root_higher,3)-centerRoot(3)))/rootYdir(2));
in_tooth = find(rootZvalue == max(rootZvalue));
point_root = rootTooth.vertices(index_root_higher(in_tooth),:);
centCrownInTooth = centerCrown + point_root - point_crown;  % �к��׼���ı�Ե����

% for debug
OBJwriteVertices('centerRoot.obj', centerRoot);  
OBJwriteVertices('point_root.obj', point_root);  
OBJwriteVertices('vz_tooth.obj', rootZvalue);  
OBJwriteVertices('��׼����y����������ĵĵ�.obj', rootTooth.vertices(index_root_higher, :));
OBJwriteVertices('�к��׼���ı�Ե����.obj', centCrownInTooth);  




%%
% 2. ���롪�����ñ�׼���Ͳ������ݵ����ĵ��Լ�������ж���
R = inv([crown_ax.x; crown_ax.y; crown_ax.z]) * axi;   % ��ת     

offset = centerPatient - centCrownInTooth * R;

tooth_T = (rootTooth.vertices) *R + repmat(offset,length(rootTooth.vertices),1);
new_crown_ax_xline = crown_ax_xline * R + repmat(offset,length(crown_ax_xline),1);
new_crown_ax_yline = crown_ax_yline * R + repmat(offset,length(crown_ax_yline),1);
new_crown_ax_zline = crown_ax_zline * R + repmat(offset,length(crown_ax_zline),1);
OBJwriteVertices('��תƽ�ƺ�ı�׼��x��.obj', new_crown_ax_xline);
OBJwriteVertices('��תƽ�ƺ�ı�׼��y��.obj', new_crown_ax_yline);
OBJwriteVertices('��תƽ�ƺ�ı�׼��z��.obj', new_crown_ax_zline);
OBJwriteVertices('��תƽ�ƺ�ı�׼������.obj', mean(tooth_T));
OBJwriteVertices('centerRoot_new.obj', centerRoot + offset);

%       2.1 ��תƽ�ƴ�����׼��
index_1 = find(patientTooth.vertex(:,2)>=(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))...
         & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - ...
         p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+0.5);
patientMergeVers = patientTooth.vertex(index_1,:);

index_2 = find(tooth_T(:,2)>=(centerPatient(2) - (patientYdir(1)*(tooth_T(:,1) - p_patient(1))+patientYdir(3)*(tooth_T(:,3) - p_patient(3)))/patientYdir(2))...
        & tooth_T(:,2)<=(p_patient(2) - (patientYdir(1)*(tooth_T(:,1) - p_patient(1))+patientYdir(3)*(tooth_T(:,3) - p_patient(3)))/patientYdir(2))+0.5);
rootMergeVers = tooth_T(index_2,:);

% for debug
 OBJwriteVertices('��һ����תƽ�ƺ�Ĵ�����׼��.obj', tooth_T);
 OBJwriteVertices('2.1�����������ں�����ĵ�.obj', patientMergeVers);
 OBJwriteVertices('2.1������׼���ں�����ĵ�.obj', rootMergeVers);

 
%       2.2 icp�㷨��һ������ 
for i = 1:length(rootMergeVers)
   [minValue,r] = mindis(patientMergeVers, rootMergeVers(i,:),1);
   minvalue(i) = minValue;  row(i) = r;
end
point11 = patientMergeVers(row,:);
[R,t,BRt,e,~,~] = icp(point11, rootMergeVers);          % icp�㷨
tooth_T = bsxfun(@plus,tooth_T*R, t);   
OBJwriteVertices('�ڶ�����תƽ�ƺ�Ĵ�����׼��.obj', tooth_T);
       


%%
% 3. �и��������ղ������ڵ��и�λ�ý����и�

%       3.1 ȷ��������׼�����и�֡�
index_root_cut =  find(tooth_T(:,2)>(p_patient(2) - ( patientYdir(1)*(tooth_T(:,1) ...
                 - p_patient(1))+ patientYdir(3)*(tooth_T(:,3) - p_patient(3)))/ patientYdir(2))-2);  %�г���������������ֵ%ע�⣡���ղ������ڷ�����
rootCutVers = tooth_T(index_root_cut ,:);


%       3.2 �ҳ����ڱ���ʱ��Ҫ�ο���ԭ�������������ϵĶ��㡣
index_patientTransform = find(patientTooth.vertex(:,2)>(centerPatient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))-1 ...
                        & patientTooth.vertex(:,2)<=(p_patient(2) - (patientYdir(1)*(patientTooth.vertex(:,1) - p_patient(1))+patientYdir(3)*(patientTooth.vertex(:,3) - p_patient(3)))/patientYdir(2))+2);
patientTransform = patientTooth.vertex(index_patientTransform,:);

% for debug
OBJwriteVertices('�ϲ�������������ֵ㼯.obj', rootCutVers);
OBJwriteVertices('�ϲ�����ο��Ĳ������ڶ���.obj', patientTransform);



%%
% 4. �ϲ�����
%       4.1 ���ں������ϲ���һ���������� 
mergedToothVers = [cutPatientVers; rootCutVers];
[mergedToothFace]=MyCrustOpen(mergedToothVers);
% for debug
writeOBJ('�ϲ�����.Obj', mergedToothVers, mergedToothFace);             % ������ ����Ƭ�������



%%
% 5. ����
mergedToothBdr = select_holes_and_boundary(mergedToothVers, mergedToothFace);
newTris = fill_mesh_holes(mergedToothVers, mergedToothFace, mergedToothBdr,'closed',200);
newTris = double(newTris);
mergeRegionIdx =  find(mergedToothVers(:,2)>(centerPatient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-1 ...
    &mergedToothVers(:,2)<=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+1);


mergeRegionVers = mergedToothVers(mergeRegionIdx,:);



newTris = reduceWrongTris(newTris);     % ȥ������Ϣ�д��������Ƭ��


writeOBJ('�����������.obj', mergedToothVers, newTris);
OBJwriteVertices('mergeRegionVers.obj', mergeRegionVers);
vecWriteUnsigned('mergeRegionIdx.obj', mergeRegionIdx);



%    �������ڲ��ֶ��һ���ֳ���������ȷ����׼�������κ���λ�ã��ҵ�tooth_T����patientTransform��������λ�ã�������Щ���Ϊ�µĵ㣬�����ڵ��λ��
indices = 1:size(mergedToothVers,1);

%    �ϲ������в���Ҫ���εĲ���
exterior = indices(mergedToothVers(:,2)<(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))-0.6...
    |mergedToothVers(:,2)>=(p_patient(2) - ( patientYdir(1)*(mergedToothVers(:,1) ...
    - p_patient(1))+ patientYdir(3)*(mergedToothVers(:,3) - p_patient(3)))/ patientYdir(2))+3);
 


%%
% 6. ����

% FOR DEBUG: ��ӡ���������
OBJwriteVertices('�ϲ�����Ķ���.obj', mergeRegionVers);
vecWriteUnsigned('�Ǻϲ�����Ķ�����������.obj', exterior);
save('�Ǻϲ�����Ķ�����������.mat', 'exterior');
OBJwriteVertices('�ϲ�����ο��Ĳ������ڶ���.obj', patientTransform);
writeOBJ('�ϲ�����', mergedToothVers, newTris);



%       6.1 �����������Ĳ���
[omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(mergedToothVers,1), newTris, exterior);
save('mergedToothVers.mat', 'mergedToothVers');
save('newTris.mat', 'newTris');
save('omega.mat', 'omega');
save('N0.mat','N0');
save('N1.mat','N1');
save('mergeRegionVers.mat','mergeRegionVers');
save('patientTransform.mat', 'patientTransform');
save('mergeRegionIdx.mat', 'mergeRegionIdx');
%[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(mergedToothVers, newTris, 'ext', 'voronoi', 'no_flatten',omega, N0, N1);

[bi_L,bi_U,bi_P,bi_Q,bi_D,bi_S,bi_M] = biharm_factor_system_modified(mergedToothVers, newTris,omega, N0, N1);



%       6.2 ִ�б��Ρ����ҵ��������ָ����ڱ��β�������ĵ㣬�������ϵĵ��������ڴ�
finalVers = mergedToothVers;
for i = 1:length(mergeRegionVers)
   [minValue,r]=mindis(patientTransform, mergeRegionVers(i,:), 1);
   minvalue(i) = minValue;  row(i) = r;
   finalVers(mergeRegionIdx(i),:) = patientTransform(r,:);
end

save('bi_L.mat', 'bi_L');
save('bi_U.mat', 'bi_U');
save('bi_P.mat', 'bi_P');
save('bi_Q.mat', 'bi_Q');
save('bi_D.mat', 'bi_D');
save('bi_S.mat', 'bi_S');
save('bi_M.mat', 'bi_M');

% finalVers = biharm_solve_with_factor( ...
%     bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
%     newTris, finalVers, omega, N0, N1, 'ext', 'no_flatten', BZ1, mergedToothVers);

finalVers = biharm_solve_with_factor_modified( ...
    bi_L, bi_U, bi_P, bi_Q, bi_D, bi_S, ...
     finalVers, omega, N0, N1);


writeOBJ('���ս��.obj', finalVers, newTris);


disp('finished');

 
