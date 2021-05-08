% dentalwithtooth1(22:28) = dentalwithtooth1(8:14);
% dentalwithtooth1(15:21) = dentalwithtooth1(8:14);
% dentalwithtooth1(8:14) =  dentalwithtooth1(1:7);
for i = 15:28
%     dentalwithtooth1(i).ID = dentalwithtooth(i).fid;
    [error,Reallignedsource,transform]=rigidICP(dentalwithtooth(i).model.vertex,...
           dentalwithtooth1(i).vertices,0,dentalwithtooth(i).model.face,dentalwithtooth1(i).faces);
%     dentalwithtooth1(i).vertices = dentalwithtooth1(i).vertices*transform.T+transform.c;
    dentalwithtooth1(i).vertices = Reallignedsource;
    figure()
    trimesh(dentalwithtooth1(i).faces,dentalwithtooth1(i).vertices(:,1),dentalwithtooth1(i).vertices(:,2),dentalwithtooth1(i).vertices(:,3))
    hold on
    trimesh(dentalwithtooth(i).model.face,dentalwithtooth(i).model.vertex(:,1),dentalwithtooth(i).model.vertex(:,2),dentalwithtooth(i).model.vertex(:,3))
end
% for n = 1:length(dentalwithtooth)/2-1
%     %%通过遍历找最近的点
%     for i = 1:length(dentalwithtooth(n).model.vertex)
%         [minValue,r]=matchest(dentalwithtooth(n+1).model.vertex,dentalwithtooth(n).model.vertex(i,:));
%         value(i) = minValue;row(i) = r;
%     end
%     k = find(value == min(value));m = row(k);
%     pointpair = [pointpair;dentalwithtooth(n).model.vertex(k,:),dentalwithtooth(n+1).model.vertex(m,:)];
%     pointfindedtooth(n,:) = mean([dentalwithtooth(n).model.vertex(k,:);dentalwithtooth(n+1).model.vertex(m,:)]); 
%     plot3(pointfindedtooth(n,1),pointfindedtooth(n,2),pointfindedtooth(n,3), 'r*','MarkerSize',100) 
%     hold on
%     value = [];
%     row = [];
% end
% for n = 1:length(dental_crown)/2-1
%     %%通过遍历找最近的点
%     for i = 1:length(dental_crown(n).model.vertex)
%         [minValue,r]=matchest(dental_crown(n+1).model.vertex,dental_crown(n).model.vertex(i,:));
%         value(i) = minValue;row(i) = r;
%     end
%     k = find(value == min(value));m = row(k);
%     pointpair = [pointpair;dental_crown(n).model.vertex(k,:),dental_crown(n+1).model.vertex(m,:)];
%     pointfindedcr(n,:) = mean([dental_crown(n).model.vertex(k,:);dental_crown(n+1).model.vertex(m,:)]); 
%     plot3(pointfindedcr(n,1),pointfindedcr(n,2),pointfindedcr(n,3), 'B*','MarkerSize',100) 
%     value = [];
%     row = [];
% end
% [d, z, tform] = procrustes(pointfindedcr, pointfindedtooth, 'Reflection',false);

for i = 1:14
    dentalwithtooth(i).model.vertex = dentalwithtooth(i).model.vertex*tform.T+tform.c(1,:);
end
for i = 15:28
    dentalwithtooth(i).model.vertex = dentalwithtooth(i).model.vertex*tform.T*[-1,0,0;0,1,0;0,0,1 ]+tform.c(1,:);
end
    
   
for i = 15:28
    fdi = dentalwithtooth1(i).ID;
    for j = 1:28
        if (upax(j).fid == fdi)
            k = j;
        end
    end
     p = mean(dentalwithtooth1(i).vertices);
%      q = mean(dental_crown(i).model.vertex);
     
     dentalwithtooth1(i).nx = upax(k).x;
     dentalwithtooth1(i).ny = upax(k).y;
     dentalwithtooth1(i).nz = upax(k).z;
     figure() 
     trimesh(dentalwithtooth1(i).faces,dentalwithtooth1(i).vertices(:,1),dentalwithtooth1(i).vertices(:,3),dentalwithtooth1(i).vertices(:,2))
%      trimesh(dentalwithtooth(i).model.face,dentalwithtooth(i).model.vertex(:,1),dentalwithtooth(i).model.vertex(:,3),dentalwithtooth(i).model.vertex(:,2))
     hold on
%      trimesh(dental_crown(i).model.face,dental_crown(i).model.vertex(:,1),dental_crown(i).model.vertex(:,3),dental_crown(i).model.vertex(:,2))
     quiver3(p(1),p(3),p(2),dentalwithtooth1(i).nx(1),dentalwithtooth1(i).nx(3),dentalwithtooth1(i).nx(2),20)
     quiver3(p(1),p(3),p(2),dentalwithtooth1(i).ny(1),dentalwithtooth1(i).ny(3),dentalwithtooth1(i).ny(2),20)
     quiver3(p(1),p(3),p(2),dentalwithtooth1(i).nz(1),dentalwithtooth1(i).nz(3),dentalwithtooth1(i).nz(2),20)
%      quiver3(q(1),q(3),q(2),nx(1),nx(3),nx(2),10)
%      quiver3(q(1),q(3),q(2),ny(1),ny(3),ny(2),10)
%      quiver3(q(1),q(3),q(2),nz(1),nz(3),nz(2),10)
end
