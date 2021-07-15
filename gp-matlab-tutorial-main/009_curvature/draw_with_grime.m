function t = draw_with_grime(vers, tris)
%Draws the input mesh with a grime effect applied to crevices of the surface.
% 模拟污垢的效果
% t = draw_with_grime(V,F);
%

%Compute the principal curvatures
kappa = discrete_curvatures(vers,tris);

%Get a function that is 0 when principal curvatures are positive, and
%increases in magnitude when they are negative
c = abs(min(0,kappa(:,1)) + min(0,kappa(:,2)));

%Plot the resulting function
t = tsurf(tris,vers, 'CData',c);
shading interp;
axis equal;
axis off;
colormap(cbrewer('Greys', 500));
light('Position',[-1.5 1 1],'Style','local');
lights = camlight;
set(t, 'FaceLighting','gouraud', 'FaceColor','interp');
set(t, 'DiffuseStrength',0.5, 'SpecularStrength',0.2, 'AmbientStrength',0.3);
camproj('perspective');
add_shadow([t],lights);

end

