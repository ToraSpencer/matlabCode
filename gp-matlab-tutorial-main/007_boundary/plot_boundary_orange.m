function plot_boundary_orange(vers, tris)
 
 
% outline()�����ҳ�����ı�Ե�ߣ�
bdryEdges = outline(tris);
bdryVerIdx = unique(bdryEdges);         %find all boundary vertices


flag = zeros(size(vers,1),1);    % �����״̬��������Ե������Ϊ1���ڲ�����Ϊ0
flag(bdryVerIdx) = 1;  


t = tsurf(tris,vers, 'CData', flag);
shading interp;
axis equal;
axis off;
cm = flipud(cbrewer('RdYlBu', 500));
colormap(cm(100:450,:));
light('Position',[-1.5 1 1],'Style','local');
lights = camlight;
set(t, 'FaceLighting','gouraud', 'FaceColor','interp');
set(t, 'DiffuseStrength',0.5, 'SpecularStrength',0.2, 'AmbientStrength',0.3);
camproj('perspective');
add_shadow([t],lights);

end

