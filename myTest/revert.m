 
%%
clear 
close all
clc

[V,F] = readOBJ('innerRevert.obj');
[OutV,OutF] = readOBJ('outer.obj');
hole = Calc_Boundary(F);
Outhole = Calc_Boundary(OutF);
%分别找出边界，然后用投影成平面问题，用triangle库补好之后，不增加点，把新增的面加入回原来的拓扑结构。
% 
sample = V(hole.boundary.edge(:,1),:);
coeff = pca(sample);

projsamp = sample*coeff;
Outprojsamp =  OutV(Outhole.boundary.edge(:,1),:)*coeff;

insize =  size(projsamp,1);
outsize = size(Outprojsamp,1);

triV = [projsamp(:,1:2);Outprojsamp(:,1:2)];

inE = [1:insize; [2:insize,1]]';
outE= [insize+1:insize+outsize; [insize+2:insize+outsize,insize+1]]';
E = [inE;outE];

H= mean(projsamp(:,1:2));

[~, newF] = triangle(triV, E, H, 'NoBoundarySteiners');



%realV = [projsamp;Outprojsamp];
%drawMesh(realV,newF);
%%
%remap
newnewF = newF;
for i=1:size(newF,1)
    for j = 1:size(newF,2)
        if newF(i,j) < 672
            newnewF(i,j) = hole.boundary.edge(newF(i,j),1);
        else
            newnewF(i,j) = size(V,1)+Outhole.boundary.edge(newF(i,j)-insize,1);
        end
    end
end

[resultV,resultF] = merge_mesh(V, F, OutV,OutF);

temp = newnewF(:,1);
newnewF(:,1) = newnewF(:,2);
newnewF(:,2) = temp;

reF = [resultF;newnewF];
drawMesh(resultV,reF);
writeOBJ('result.obj',resultV,reF);
    view(3)
    axis equal
    axis off
    camlight
    lighting gouraud
    cameratoolbar
    set(gca, 'Position',[0 0 1 1]);
 