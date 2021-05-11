% clc
% clear all
% 
% all = Read_Obj('output_roots.obj');
vert1 = [];vert2 = [];vert3 = [];
n_ya = 14;
houzhui = '.stl';
figure()
for num =5:n_ya
    serFile = num2str(num);
    filename = strcat(serFile,houzhui);
    fv = stlread(filename);
    v_sym = [-fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3)];
    vertices = [fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3)];
    vert1 = [vert1;vertices(fv.faces(:,1),:)];
    vert2 = [vert2;vertices(fv.faces(:,2),:)];
    vert3 = [vert3;vertices(fv.faces(:,3),:)];
%     trimesh(fv.faces,vertices(:,1),vertices(:,2),vertices(:,3))
    [~,MC]=curvatures(vertices(:,1),vertices(:,2),vertices(:,3),fv.faces);
%     hold on 
    pp= patch('Faces',fv.faces,'Vertices',vertices,'FaceVertexCData',MC,'FaceColor','interp','EdgeColor','none');
     caxis([-25 , 24])
     colormap jet
     colorbar
end

% figure()
% pointpair = [];
% for n = 6:33
%      face = [];
%      l = length(all(n,1).vertex);
%      face = [l+all(n,1).face(:,3)+1,l+all(n,1).face(:,2)+1,l+all(n,1).face(:,1)+1];
% %      trimesh(face,all(n,1).vertex(:,1),all(n,1).vertex(:,2),all(n,1).vertex(:,3))
%      [~,MC]=curvatures(all(n,1).vertex(:,1),all(n,1).vertex(:,2),all(n,1).vertex(:,3),face);
%      hold on
%      pp= patch('Faces',face,'Vertices',all(n,1).vertex,'FaceVertexCData',MC,'FaceColor','interp','EdgeColor','none');
%      caxis([-4 , 4])
%      colormap jet
%      colorbar
%    
% end