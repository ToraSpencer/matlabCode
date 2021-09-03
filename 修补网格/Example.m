clc
clear all
close all

% load date
load('sphere_8k.mat');

src = 1;

D = perform_fast_marching_mesh(surfdata.point, surfdata.face, src);

% 球上待求顶点索引
idx_unknow = find(D <= 1);

% 球上已知顶点索引
idx_know = find(~ismember(1:size(surfdata.point,1), idx_unknow));

A = adjacency_matrix(surfdata.face);
L = speye(size(A)) - spdiags(1./sum(A,2), 0, size(A,1), size(A,1)) * A;

V = solve_equation(L, zeros(size(surfdata.point)), idx_know, surfdata.point(idx_know,:));
F = surfdata.face;

% 待求问题图示
figure
drawMesh(V, F, 'facecolor','y', 'edgecolor','none');
plot3(V(idx_know,1), V(idx_know,2), V(idx_know,3), 'r.')
plot3(V(idx_unknow,1), V(idx_unknow,2), V(idx_unknow,3), 'b.')

view(3);
axis equal
axis off
camlight
lighting gouraud
set(gca, 'Position',[0 0 1 1]);

% 优化求解
L = cotmatrix(V, F);
M = massmatrix(V, F);

A = L * inv(M) * L;% * inv(M) * L;
%A=A*A*A;
A = 0.5*blkdiag(A, A, A);
B = zeros(3*size(V,1),1);
Aeq = sparse(1:3*length(idx_know), [idx_know,idx_know+size(V,1),idx_know+2*size(V,1)], 1, 3*length(idx_know), 3*size(V,1));
Beq = reshape(V(idx_know,:), [], 1);

X = quad_prog(A, B, Aeq, Beq, zeros(0,3*size(V,1)), zeros(0,1));

Vnew = reshape(X, [], 3);

figure
drawMesh(Vnew, F, 'facecolor','y', 'edgecolor','none');
plot3(Vnew(idx_know,1), Vnew(idx_know,2), Vnew(idx_know,3), 'r.')
plot3(Vnew(idx_unknow,1), Vnew(idx_unknow,2), Vnew(idx_unknow,3), 'b.')

view(3);
axis equal
axis off
camlight
lighting gouraud
set(gca, 'Position',[0 0 1 1]);