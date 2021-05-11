% clc
% clear all
% % houzhui = '.stl';
functionname='test.m'; functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir 'trianglerayintersection'])
addpath([functiondir 'toothmodel'])
addpath([functiondir 'theratofp'])
addpath([functiondir 'Root deformation'])
addpath([functiondir 'modelread'])
addpath([functiondir 'mixFE'])
addpath([functiondir 'margin line'])
addpath([functiondir 'Arch curve'])
addpath([functiondir '��׼ģ��'])

% vert1 = [];vert2 = [];vert3 = [];
% % center1 = zeros(n_ya,3);center2 = zeros(n_ya,3);
all = Read_Obj('117608_lowerMesh.obj');
n_ya = length(all);
%%1.��ȡ���˵�����
yazhou = textread('117608_lowerDir_centre.txt');
for i = 1:n_ya
    case_center(i,:) = yazhou((i-1)*4+1,:);
    case_axle{i} = yazhou((i-1)*4+2:i*4,:);
end
%%2.��ȡ���ݵı��

% load A
% figure()
pointpair = [];
for n = 1:n_ya-1
     face = [];
     l = length(all(n,1).vertex);
     face = [l+all(n,1).face(:,3)+1,l+all(n,1).face(:,2)+1,l+all(n,1).face(:,1)+1];
       
%      %����͵�
     num = find(all(n,1).vertex(:,3) == min(all(n,1).vertex(:,3)));
     Hoz(n,:) = all(n,1).vertex(num,:);  
%     trimesh(face,all(n,1).vertex(:,1),all(n,1).vertex(:,2),all(n,1).vertex(:,3))
%      hold on
    %%ͨ������������ĵ�
    for i = 1:length(all(n,1).vertex)
        [minValue,r]=matchest(all(n+1,1).vertex,all(n,1).vertex(i,:));
        value(i) = minValue;row(i) = r;
    end
    k = find(value == min(value));m = row(k);
    pointpair = [pointpair;all(n,1).vertex(k,:),all(n+1,1).vertex(m,:)];
    pointfinded(n,:) = mean([all(n,1).vertex(k,:);all(n+1,1).vertex(m,:)]); 
%     plot3(pointfinded(n,1),pointfinded(n,2),pointfinded(n,3), 'r*','MarkerSize',100) 
    value = [];
    row = []; 
   
end
%%��xȫ��ͼ�ϲ�����
I  = imread('12.jpg');
figure()
imshow(I)
[x,y] = ginput(length(pointfinded));
close;
mid = floor(length(pointfinded)/2);
for i =1:length(pointfinded)-1
   dis1(i) = sqrt((x(i)-x(i+1))^2+(y(i)-y(i+1))^2);
   dis(i) = sqrt(sum((pointfinded(i,:)-pointfinded(i+1,:)).^2));    
end


rat = sum(dis)/sum(dis1);
Hy = y(mid)-y;Hz = pointfinded(:,3)- pointfinded(mid,3);xlong = abs(pointfinded(mid,1)- max(pointfinded(:,1)));
Hy = Hy*rat;
thetaend = asin(Hy/xlong);thetastar = asin(Hz/xlong);
%%����������
vert1 = [];vert2 = [];vert3 = [];
figure()
for n = 1:33
    face = [];
     l = length(all(n,1).vertex);
     face = [l+all(n,1).face(:,3)+1,l+all(n,1).face(:,2)+1,l+all(n,1).face(:,1)+1];
     Face{n} = face;
     [ver] = rotation_xyz1(all(n,1).vertex,0,0.4,0,Hoz(mid,:),0);
     v_sym = [ver(:,2),ver(:,1),ver(:,3)]; 
     vert1 = [vert1;v_sym(face(:,1),:)];
     vert2 = [vert2;v_sym(face(:,2),:)];
     vert3 = [vert3;v_sym(face(:,3),:)];
    
     center(n,:) = mean(v_sym);
     trimesh(face,v_sym(:,1),v_sym(:,2),v_sym(:,3))
     hold on
     Nor{n} = per_vertex_normals(all(n,1).vertex,face); 
end 
[A,N]= pfit(center);
N = 3;
A = polyfit(center(:,1),center(:,2),N); 
%���������
x0=40:-0.1:-40;
y0=A(1)*x0.^N    ; %������õ�ϵ����ʼ�����������ʽ����
for num=2:1:N+1     
    y0=y0+A(num)*x0.^(N+1-num);    
end
plot(x0,y0)
% for x0 =-30:2:30
% y0 = A(1)*x0.^3+A(2)*x0.^2+A(3)*x0+A(4);
% end
NN = 400;%�����õ���
pt =  interparc(NN,x0,y0,'spline');
FX = 3*A(1)*pt(:,1).^2+2*A(2)*pt(:,1)+A(3);
FY = -ones(NN,1);










%% �ڹ��߷����ϵ���������
%1.ȷ����άģ�͸�����λ��
%2.�������ȫ��ͼ�ϵ�ͶӰ�����е�ȫ��ͼ�Ա�
%%2.1�������������Ӵ����������ĵ㣬.
%%2.2�ҵ�����㣬��άģ���ϵ�
%3.ȷ��������λ�ã��ص�������ô�������ȫ��ͼ��ʵ�ʵ�ȫ��ͼ���Ա�
%���⣺���ģ�͵�������ʵ�ʵ������Σ�Ӧ����Σ�
%�漰��������ô��ȫ��ͼ��ͼƬͶӰ���ο�ƽ����

% figure()
% imshow(I)
% [x,y] = ginput(6);
% points = {[x(1) y(1)],[x(2) y(2)],...
%            [x(3) y(3)],[x(4) y(4)],...
%            [x(5) y(5)],[x(6) y(6)],...
%            };
%  hobbysplines(points,'debug',true,'cycle',false);
for n = 28:33
    %��׼ģ�ͺͲ����õ������ڣ��ñ�ģ���ղ��˵����ڸı䷽��
    
    
    
    
    
    
    
    
    
    
    
    
    
    j =1;
    face = [];
    face = Face{1,n};
%     figure()
%     trimesh(face,all(n,1).vertex(:,1),all(n,1).vertex(:,2),all(n,1).vertex(:,3))
%     hold on
    %Ѱ�Ҹ���㣺��͵�p_apx
    n_apx = find(all(n,1).vertex(:,3)==min(all(n,1).vertex(:,3)));
    p_apx = all(n,1).vertex(n_apx,:);
    %Ѱ����ά���������ĵ�thcenter
    Po = all(n,1).vertex(find(all(n,1).vertex(:,3)<=-12 &all(n,1).vertex(:,3)>=-13),:);
    thcenter = mean(Po);
    %��Ӧ��ȫ��ͼ�ϵ�����
    figure()
    imshow(I)
    %��ȡ�����Ϊ�����㣬���м��Ϊ������
    Npick = 13;
    [x,y] = ginput(Npick);
    close
    %����ÿ�����Ӧ��ʵ������
    %��ά�ϵ����ĵ�
    twcenter = [mean([x(1),x(end)]) mean([y(1),y(end)])];
    %ÿ����������ĵ��λ������
    twpoint_center = [x y] - twcenter*rat;
    %��ά�϶�Ӧ��������λ��,���ǽ���ά�����ƽ���ֳɣ�Npick-1��/2-1��
    %��z�ĸ߶�ƽ������
    deltaz = (p_apx(3) - thcenter(3))/((Npick-1)/2-1);
    V = all(n,1).vertex;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
    rest_V = V;
% %     %Ѱ�ҿ��Ʊ仯�ĵ㣬�������仯λ��
%     %1.Ѱ�����д�ֱ�����ߵĵ㣬����������ƽ���ϵĵ�
%       for i = 1:length(rest_V)
%         d = sqrt((pt(:,1)-rest_V(i,2)).^2 + (pt(:,2)-rest_V(i,1)).^2); %����С���벻�Ǻܺ��ʣ����ĳ�����ݲ��������������ؾ� 
%         n_mindis = find(d == min(d));
%         r = [FX(n_mindis),FY(n_mindis),0];
%         if abs(dot(r,Nor{1,n}(i,:)))<=0.04
%             Edge_point{n}(j,:) = [rest_V(i,:)];
%             dotvalue{n}(j) = abs(dot(r,Nor{1,n}(i,:)));
%             j =j+1;
%         end       
%      end
%     Edge_points = [];pf = Edge_point{n};Edge_points = pf(find(pf(:,3)<=thcenter(3)),:);
%     %2.Edge_points(:,3)����ÿ�θ߶���Ϊ�ӽ��ĵ㷢��λ�ñ仯
    deltaH = (p_apx(3)-thcenter(3))/((Npick-1)/2-1);
bi_V = zeros(size(V));
re_V = rest_V;
bi_bndtype = 'ext';
BZ1 = zeros(size(V,1),3);
reduction = 'no_flatten';masstype = 'voronoi';
indices = 1:size(V,1);
exterior = indices( V(:,3) > -12 | V(:,3) <= min(V(:,3))+0.1);
[Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(V,1), face, exterior);

[bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(V,face, bi_bndtype,masstype,reduction,Omega,N0,N1);
%��z��Ĵ�Сֵ

V_changed  = z_change(re_V, y,thcenter,p_apx,rat,1);
% figure()
    for n_z = 1:(Npick-1)/2-1
        n_point = [];
        HZ(n_z) = thcenter(3)+(n_z-1)*deltaH;
% %         H_dis = abs(Edge_points(:,3)-HZ);
% %         n_Edegepoint(n_z) = find(H_dis == min(H_dis));
% %         control_point{n_z} = Edge_points(n_Edegepoint(n_z),:);
        n_point = find(re_V(:,3)>=HZ(n_z)-1 & re_V(:,3)<=HZ(n_z)+1);
        circle_point{n_z} = re_V(n_point,:);
        centerofcircle = mean(circle_point{n_z});
        d = sqrt((pt(:,1)-centerofcircle(2)).^2 + (pt(:,2)-centerofcircle(1)).^2);
        n_dmin = find(d == min(d));
        if pt(n_dmin,2) < centerofcircle(1)
            ints = [pt(:,1)-min(d).*FX./(FX.^2+FY.^2).^0.5,pt(:,2)-min(d).*FY./(FX.^2+FY.^2).^0.5];%ints:��������Ǩ�Ƶ��������ĵ㣬������Ϊ���Ĳο�
        else
            ints = [pt(:,1)+min(d).*FX./(FX.^2+FY.^2).^0.5,pt(:,2)+min(d).*FY./(FX.^2+FY.^2).^0.5];
        end
        figure()
%         plot(ints(:,2),ints(:,1))
%         hold on 
%         plot(centerofcircle(1),centerofcircle(2),'*r')
%         plot(circle_point{n_z}(:,1),circle_point{n_z}(:,2))
%         hold on
        %Ѱ�Ҹ����߽ӽ�ƽ�еĵ�
        plot3(circle_point{n_z}(:,1),circle_point{n_z}(:,2),circle_point{n_z}(:,3),'*g')
        hold on
        for k = 1:length(n_point)
            Ndir  = circle_point{n_z}(k,1:2) - [ints(:,2),ints(:,1)];%�ƶ������������ϵĵ���ĳһ�����Ȧ��֮�������
            dotvalue = sum([FY,FX].* Ndir,2);
            theta_ray = acos(dotvalue./(sqrt(FX.^2+FY.^2).*sqrt(Ndir(:,1).^2+Ndir(:,2).^2)));%�н�
            angle_min_loc(k) = find(theta_ray == min(theta_ray));%������ƽ�е�����
            angle_min(k) = min(theta_ray);        %��С�ļн�ֵ       
        end
        n_dot = find(ints(angle_min_loc,1) <= max(circle_point{n_z}(:,2))&ints(angle_min_loc,1) >= min(circle_point{n_z}(:,2)));%�������߹���y��Գ�
        
        dot_circle = circle_point{n_z}(n_dot,:);
%         plot3(dot_circle(:,1),dot_circle(:,2),dot_circle(:,3),'*r')
%         plot(ints(:,2),ints(:,1))
        margdotloc = angle_min_loc(n_dot);
        margdot = ints(margdotloc,1:2);
        dot_circle_max = dot_circle(find(dot_circle(:,1)>= centerofcircle(1)),:);dot_circle_min = dot_circle(find(dot_circle(:,1)< centerofcircle(1)),:);
        dotmax = [];
        dotmax = margdot(find(margdotloc>=n_dmin),:);dotmin = margdot(find(margdotloc<n_dmin),:);
        dmax = [];dmin = [];
%         for i = 1:length(dotmax)
%             dmax(i) =  sqrt((dot_circle_max(i,1)-dotmax(i,2)).^2 + (dot_circle_max(i,2)-dotmax(i,1)).^2);      
%         end
%          for i =1:length(dotmin)
%             dmin(i) =  sqrt((dot_circle_min(i,1)-dotmin(i,2)).^2 + (dot_circle_min(i,2)-min(i,1)).^2);      
%         end
%         ndmax = find(abs(dmax) == min(abs(dmax)));
%         ndmin = find(abs(dmin) == min(abs(dmin)));
        [ndmax,~] = mindispoint(dot_circle_max(:,1:2),ints);
        [ndmim,~] = mindispoint(dot_circle_min(:,1:2),ints);
        base = [dot_circle_max(ndmax,:);dot_circle_min(ndmin,:)];      
%         plot(margdot(:,2),margdot(:,1),'*b')      
        shape = [];
%         plot(margmax(:,1),margmax(:,2),'*k')
%         plot(margmin(:,1),margmin(:,2),'*m')
        for i = n_dmin:length(pt)
            dfidmax(i-n_dmin+1) =  sqrt((thcenter(1)-ints(i,2)).^2 + (thcenter(2)-ints(i,1)).^2) - abs((x(n_z)-(x(1)+x(Npick))/2)*rat);      
        end
         for i =1:n_dmin
            dfidmin(i) =  sqrt((thcenter(1)-ints(i,2)).^2 + (thcenter(2)-ints(i,1)).^2) - abs((x(Npick+1-n_z)-(x(1)+x(Npick))/2)*rat);      
        end
        ndfidmax = find(abs(dfidmax) == min(abs(dfidmax)));
        ndfidmin = find(abs(dfidmin) == min(abs(dfidmin)));  
% %         ���ֶ����������������ѡȡ�������������ĵ�  
        shape(1,:) = [ints(n_dmin+ndfidmax-1,2),ints(n_dmin+ndfidmax-1,1),(y(n_z)-(y(1)+y(Npick))/2)*rat+centerofcircle(3)];  
        shape(2,:) = [ints(ndfidmin,2),ints(ndfidmin,1),(y(Npick+1-n_z)-(y(1)+y(Npick))/2)*rat+centerofcircle(3)];
        plot(shape(:,1),shape(:,2),'*k')
        plot(base(:,1),base(:,2),'*m')
        %%�õ�������ƣ��ᵥ���γ�С͹�⣬���Ǹı�һȦ���λ�ã�����ʹ�����λ��
        %1.�������λ��ȷ��ƽ�ƺ����ű�׼
        [~, z, tform] = procrustes( shape(:,1:2),base(:,1:2));%�˴�baseΪ������ԭʼ�ĵ�
       %% Z = b*Y*T + c;
        tranpoints =tform.b * circle_point{n_z}(:,1:2) * tform.T + tform.c(1,:);
        plot3(tranpoints(:,1),tranpoints(:,2),circle_point{n_z}(:,3),'*c')
        hold on
        V(n_point,1:2) = tranpoints;
%         BZ1 = zeros(size(V,1),3);
%         if(strcmp(bi_bndtype,'deriv'))
%           % grab just inner ring boundary
%           inner_N0 = N0( (V(N0,1)-0.5).^2 + (V(N0,2)-0.5).^2 <= 0.2^2);
%           % grad just outer ring boundary
%           outer_N0 = N0( (V(N0,1)-0.5).^2 + (V(N0,2)-0.5).^2 >= 0.4^2);
% 
%           % magnitude of tangents, roughly propotional to edge length
%           inner_magnitude = 0.1;
%           outer_magnitude = 0.05;
% 
%           % Inner ring
%           % point tangents away from center
%           BZ1(inner_N0,:) = ...
%             [(V(inner_N0,1)-0.5) (V(inner_N0,2)-0.5) 0.0*V(inner_N0,3)];
%           % normalize and multiply by magnitude
%           BZ1(inner_N0,:) = ...
%             BZ1(inner_N0,:)./repmat(sqrt(sum(BZ1(inner_N0,:).^2,2)),1,3)*inner_magnitude;
% 
%           % Outer ring
%           % point tangents toward from center
%           BZ1(outer_N0,:) = ...
%             [(0.5-V(outer_N0,1)) (0.5-V(outer_N0,2)) 0.0*V(outer_N0,3)];
%           % normalize and multiply by magnitude
%           BZ1(outer_N0,:) = ...
%             BZ1(outer_N0,:)./repmat(sqrt(sum(BZ1(outer_N0,:).^2,2)),1,3)*outer_magnitude;
%         end
%         bi_V = biharm_solve_with_factor( bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, face, V, Omega, N0, N1, bi_bndtype, reduction,BZ1,re_V);
%         re_V = bi_V;
        
        clear margdotloc
        clear angle_min
        clear n_dot
        clear margdot0fcircle
        clear n_point
        clear angle_min
        clear angle_min_loc
        clear dfidmax
        clear dfidmin
        

        
        
        
    end    
    
    
    
    
    
    
%     %������������������
    l_pro = sqrt((x(4)-mean([x(1),x(end)]))^2+(y(4)-mean([y(1),y(end)]))^2)*rat;%�˴����ݺ�����Ƿ�Ҫ�ֿ���
%     %��άģ������������
    l = sqrt(sum((p_apx-thcenter).^2));
    %�ָ��ⳤ�ϳ����߸��ⳤ�Ƚ϶̣����˾��ò�Ҫ������������Ǳ任
    %��ʱ���Բ��ÿ����õ����ķ���
    %����ʹ�û������Ԫ�ķ���
    %1.V��ȷ��,���������еĵ�
   
    bi_V = zeros(size(V));
    bi_bndtype = 'ext';
    BZ1 = zeros(size(V,1),3);
    reduction = 'no_flatten';masstype = 'voronoi';
    indices = 1:size(V,1);
    exterior = indices( V(:,3) > -12);
    [Omega, N0, N1, N2, outside_region_of_interest] = layers_from_handle(size(V,1), face, exterior);
   
    [bi_L,bi_U,bi_P,bi_Q,bi_R,bi_S,bi_M] = biharm_factor_system(V,face, bi_bndtype,masstype,reduction,Omega,N0,N1);
    
    %2.������������λ�ñ仯
    %�����������ķ���������������᷽�����˶�,�����������λ��
    r = p_apx - thcenter;
    r = r/sqrt(sum(r.*r));
    
    V(n_apx,:) = [l_pro*(p_apx(1)-thcenter(1))/l+thcenter(1),...
        l_pro*(p_apx(2)-thcenter(2))/l+thcenter(2),...
        l_pro*(p_apx(3)-thcenter(3))/l+thcenter(3)];
  
    bi_V = biharm_solve_with_factor( ...
        bi_L, bi_U, bi_P, bi_Q, bi_R, bi_S, bi_M, ...
        face, V, Omega, N0, N1, bi_bndtype, reduction,BZ1,rest_V);
    
    %1.ȷ���������Ϲ̶�λ��uf����Ե��
    
    
    
   
% %     if l>l_pro
%         %��ͶӰ�����ϵı仯
%         %������ת�ĽǶ�
%         theta_R = 0.2;%angle(acos ((abs(l)-abs(abs(p_apx(3)-p_cen(3))-abs(y(4)-mean([y(1),y(end)]))))/abs(l)));
%         %ͨ���������Ĺ��ߵķ���,����ͨ������������������ϵ�����������
%         s = sqrt((p_apx(1)-pt(:,1)).^2+(p_apx(2)-pt(:,2)).^2);
%         g = find (s == min(s));
%         % ��������ת������cenΪ���ģ���lcen-apx��ray��ɵ�ƽ������ת��
%         lcen_apx = p_apx-thcenter; 
%         n_ray = [p_apx(1)-pt(g,1),p_apx(2)-pt(g,2),0];
%         dirct = cross(lcen_apx,n_ray);
%         kk = find(all(n,1).vertex(:,3)<=thcenter(3));
%         P = all(n,1).vertex(kk,:);
%         Pr=rot3d1(P,thcenter,dirct,theta_R);
%         PP = all(n,1).vertex;
%         for i = 1:length(kk)
%             PP(kk(i),:) = Pr(i,:);
%         end
%         figure()
%         h=trisurf(face,PP(:,1),PP(:,2),PP(:,3),0.3,'Edgecolor','none');
%         hold
%         light
%         lighting phong;
%         set(gca, 'visible', 'off')
%         set(gcf,'Color',[1 1 0.88])
%          h=trisurf(face,all(n,1).vertex(:,1),all(n,1).vertex(:,2),all(n,1).vertex(:,3),0.8,'Edgecolor','none');
%          
%     else
        
%     end
%     %%��˹ţ�ٵ���
%     E =  sum ((P+d-Pr).^2);
    
    
    
    
%     a = 
    
    
    
    for i = 1:length(all(n,1).vertex)
        d = sqrt((pt(:,1)-all(n,1).vertex(i,2)).^2 + (pt(:,2)-all(n,1).vertex(i,1)).^2); 
        n_mindis = find(d == min(d));
        r = [FX(n_mindis),FY(n_mindis),0];
%         [n_clo,n_nu] = find(face == i);
%         n_loc = max(n_clo);
        if abs(dot(r,Nor{1,n}(i,:)))<=0.05
            Edge_point{n}(j,:) = [all(n,1).vertex(i,:)];
            dotvalue{n}(j) = abs(dot(r,Nor{1,n}(i,:)));
            
            plot3(all(n,1).vertex(i,1),all(n,1).vertex(i,2),all(n,1).vertex(i,3),'*r')
            hold on
            j =j+1;
        end       
    end
    Edge_points = [];pf = Edge_point{n};Edge_points = pf(find(pf(:,3)<=-12.53),:);%�˴���Ҫ��Ե��Ϣ
    n = fitNormal(Edge_points,1);
 
    hobbysplines(points,'debug',true,'cycle',false);
    view(3)
    hold on
    clear points
end
       