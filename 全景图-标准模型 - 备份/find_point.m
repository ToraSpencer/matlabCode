clc
clear all
n_ya = 14;
% vert1 = [];vert2 = [];vert3 = [];
addpath = 'C:\Users\lishuping\Desktop\CT仿真\TriangleRayIntersection\TriangleRayIntersection\未修补牙齿网格(1)\未修补牙齿网格\l_1';
% center1 = zeros(n_ya,3);center2 = zeros(n_ya,3);
% all_struct = Read_Obj('output_roots.obj');

figure()
% nextfilePath = '.obj';
for count = 0:n_ya-1
    serFile = num2str(count);
    fileName = ['splittooth_',serFile,'.obj'];
%     fullfilePath = strcat(fileName);
    [preVerArray ,preFacArray] =  ReadObj(fileName);
%     vert1 = [vert1;preVerArray(preFacArray(:,1),:)];
%     vert2 = [vert2;preVerArray(preFacArray(:,2),:)];
%     vert3 = [vert3;preVerArray(preFacArray(:,3),:)];
    all(count+1,1).preVerArray = preVerArray;
    all(count+1,1).vertex1 = preVerArray(preFacArray(:,1),:);
    all(count+1,1).vertex2 = preVerArray(preFacArray(:,2),:);
    all(count+1,1).vertex3 = preVerArray(preFacArray(:,3),:);
    all(count+1,1).face = preFacArray;
    all(count+1,1).center =mean(preVerArray);
    trimesh(preFacArray,preVerArray(:,1),preVerArray(:,2),preVerArray(:,3))
    hold on 
end

pointpair = [];
for n = 1:n_ya-1
    %%通过遍历找最近的点
    for i = 1:length(all(n,1).preVerArray)
        [minValue,r]=matchest(all(n+1,1).preVerArray,all(n,1).preVerArray(i,:));
        value(i) = minValue;row(i) = r;
    end
    k = find(value == min(value));m = row(k);
    pointpair = [pointpair;all(n,1).preVerArray(k,:),all(n+1,1).preVerArray(m,:)];
    pointfinded(n,:) = mean([all(n,1).preVerArray(k,:);all(n+1,1).preVerArray(m,:)]); 
    plot3(pointfinded(n,1),pointfinded(n,2),pointfinded(n,3), 'r*','MarkerSize',100)
%     plot3(all(n,1).preVerArray(k,1),all(n,1).preVerArray(k,2),all(n,1).preVerArray(k,3), 'r*','MarkerSize',100)
%     plot3(all(n+1,1).preVerArray(m,1),all(n+1,1).preVerArray(m,2),all(n+1,1).preVerArray(m,3),'r*','MarkerSize', 100)
    value = [];
    row = [];
        
    
    %%通过重心找与牙齿的相交点，磨牙部分会出现角点在磨牙的牙凹处
%     plot3(all(n+1,1).center(:,1),all(n+1,1).center(:,2),all(n+1,1).center(:,3),'*r')
%     plot3(all(n,1).center(:,1),all(n,1).center(:,2),all(n,1).center(:,3),'*r')
%     dir = all(n+1,1).center-all(n,1).center;
%     orig = all(n,1).center;
%     c = sqrt(sum(dir.^2));
%     [intersect, t, u, v, xcoor] = TriangleRayIntersection(orig, dir, [all(n,1).vertex1;all(n+1,1).vertex1], ...
%                                   [all(n,1).vertex2;all(n+1,1).vertex2], [all(n,1).vertex3;all(n+1,1).vertex3],...
%                               'lineType','segment','planeType','one sided');
%      scatter3(xcoor(intersect,1), xcoor(intersect,2), xcoor(intersect,3), 100, 'b', 'o', 'filled')
%      line('XData',orig(1)+[0 dir(1)],'YData',orig(2)+[0 dir(2)],'ZData',...
%        orig(3)+[0 dir(3)],'Color','r','LineWidth',1)
%       


end
% n_ya = 14;%length(all_struct);   
% for n =1:n_ya
%     face = [];
%     l = length(all_struct(n,1).vertex);
%      face = [l+all_struct(n,1).face(:,3)+1,l+all_struct(n,1).face(:,2)+1,l+all_struct(n,1).face(:,1)+1];
%     vertices = [all_struct(n,1).vertex(:,2),all_struct(n,1).vertex(:,1),all_struct(n,1).vertex(:,3)];
% %     vertices = [fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3)];
%     vert1 = [vert1;vertices(face(:,1),:)];
%     vert2 = [vert2;vertices(face(:,2),:)];
%     vert3 = [vert3;vertices(face(:,3),:)];
%     trimesh(face,vertices(:,1),vertices(:,2),vertices(:,3))
%    
%     hold on 
% end
j =1;
c = 45;
dir = [0,1,0];
for x = -40:0.2:40
    i =1;
    for z = -40:0.2:40
        orig = [x,-10,z];
        [intersect, t, u, v, xcoor] = TriangleRayIntersection(orig, dir, vert1, vert2, vert3);
        T = find(~isnan(xcoor(:,2)));
        if isempty(T)~=1
            Spix(j,i) = dir(2)*(T(end)-orig(2)); 
        else
            Spix(j,i)=0;
        end
%         if sum(intersect)~=0
%             d = t(find(intersect~=0));
%             if mod(sum(intersect),2)==0
%                 
%                 Spix(j,i) = abs(d(end)); 
%             end
%         else
%             Spix(j,i)=0;
%         end
        i = i+1;
%          plot3(orig(1),orig(2),orig(3),'*r')
%          scatter3(xcoor(intersect,1), xcoor(intersect,2), xcoor(intersect,3), 100, 'b', 'o', 'filled')
%      line('XData',orig(1)+[0 dir(1)*c],'YData',orig(2)+[0 dir(2)*c],'ZData',...
%          orig(3)+[0 dir(3)*c],'Color','r','LineWidth',1)
    end
%     fprintf('intresections=%i; time=%f sec\n', sum(intersect), toc) 
    
    j = j+1;
end
figure()
imshow(Spix/max(max(Spix)))
% Spix = Spix*256/max(max(Spix));
% [Fx,Fy] = gradient(Spix);
% Theta = atan(Fy./Fx);
% mx = Fx./sqrt(Fx.^2+Fy.^2);
% my = Fy./sqrt(Fx.^2+Fy.^2);
% [sx,sy] = gradient(my);
% [cx,cy] = gradient(mx);
% D = sqrt(sx.^2+sy.^2+cx.^2+cy.^2);
% lapla = calDifLaplacian(Spix,'standard');
% D1 = zeros(size(Spix));D2 = zeros(size(Spix));
% D1(find(lapla<0)) = D(find(lapla<0));
% D2(find(Spix>0)) = D1(find(Spix>0));
% t = mean(D2(find(D2>0)))+std(D2(find(D2>0)));
% D2(find(D2<t)) = 0;
% [m,n] = find(D2>t);

