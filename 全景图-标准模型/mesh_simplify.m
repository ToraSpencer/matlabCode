%% surface simplification using quadratic error metrics
% clc
% clear all

% 
% load('all.mat');

% V = surfdata.point;
% F = surfdata.face;
A = all(1,:);
V = A.vertex(:,1:3);
l = length(A.vertex);
F = [l+A.face(:,3)+1,l+A.face(:,2)+1,l+A.face(:,1)+1];
nv = size(V,1); % total vertex number
np = 0.1*nv; % remained vertex number

% fundamental error quadric
N = normals(V,F);
N = normalizeVector3d(N);
p = [N, -sum(N .* V(F(:,1),:), 2)];
Q0 = bsxfun(@times, permute(p, [2,3,1]), permute(p, [3,2,1]));

% compute the Q matrices for all the initial vertices.
nf = size(F,1);
Q = zeros(4,4,nv);
valence = zeros(nv,1);
for i = 1:nf
    for j = 1:3
        valence(F(i,j)) = valence(F(i,j)) + 1;
        Q(:,:,F(i,j)) = Q(:,:,F(i,j)) + Q0(:,:,i);
    end
end

% all valid pairs
TR = triangulation(F,V);
E = edges(TR);

% compute Q1+Q2 for each pair
Qbar = Q(:,:,E(:,1)) + Q(:,:,E(:,2));

% a simple scheme: select either v1, v2 or (v1+v2)/2
ne = size(E,1);
v1 = permute([V(E(:,1),:),ones(ne,1)], [2,3,1]);
v2 = permute([V(E(:,2),:),ones(ne,1)], [2,3,1]);
vm = 0.5 .* (v1 + v2);
v = [v1, v2, vm];

cost = zeros(ne,3);
cost(:,1) = sum(squeeze(sum(bsxfun(@times,v1,Qbar),1)).*squeeze(v1),1)';
cost(:,2) = sum(squeeze(sum(bsxfun(@times,v2,Qbar),1)).*squeeze(v2),1)';
cost(:,3) = sum(squeeze(sum(bsxfun(@times,vm,Qbar),1)).*squeeze(vm),1)';

figure;
draw_t = 0;
num = nv;
tic
for i = 1:nv-np
    if (nv - i) < 0.9*num
        num = nv - i;
        
        clf
        drawMesh(V, F, 'facecolor','y', 'edgecolor','k', 'linewidth', 1.2);
        view([0 90])
        axis equal
        axis off
        camlight
        lighting gouraud
        cameramenu
        drawnow
    end
    
    [min_cost, vidx] = min(cost,[],2);
    [useless, k] = min(min_cost);
    e = E(k,:);

    % update position for v1
    V(e(1),:) = v(1:3, vidx(k), k)';
    V(e(2),:) = nan;

    % update Q for v1
    Q(:,:,e(1)) = Q(:,:,e(1)) + Q(:,:,e(2));
    Q(:,:,e(2)) = nan;

    % updata face
    F(F == e(2)) = e(1);
    f_remove = sum(diff(sort(F,2),[],2) == 0, 2) > 0;
    F(f_remove,:) = [];

    % collapse and delete edge and related edge information
    E(E == e(2)) = e(1);
    E(k,:) = [];
    cost(k,:) = [];
    Qbar(:,:,k) = [];
    v(:,:,k) = [];

    % delete duplicate edge and related edge information
    [E,ia,ic] = unique(sort(E,2), 'rows');
    cost = cost(ia,:);
    Qbar = Qbar(:,:,ia);
    v = v(:,:,ia);
 
    % pairs involving v1
    pair = sum(E == e(1), 2) > 0;
    npair = sum(pair);

    % updata edge information
    Qbar(:,:,pair) = Q(:,:,E(pair,1)) + Q(:,:,E(pair,2));
    
    pair_v1 = permute([V(E(pair,1),:),ones(npair,1)], [2,3,1]);
    pair_v2 = permute([V(E(pair,2),:),ones(npair,1)], [2,3,1]);
    pair_vm = 0.5 .* (pair_v1 + pair_v2);
    v(:,:,pair) = [pair_v1, pair_v2, pair_vm];
    
    cost(pair,1) = sum(squeeze(sum(bsxfun(@times,pair_v1,Qbar(:,:,pair)),1)).*squeeze(pair_v1),1)';
    cost(pair,2) = sum(squeeze(sum(bsxfun(@times,pair_v2,Qbar(:,:,pair)),1)).*squeeze(pair_v2),1)';
    cost(pair,3) = sum(squeeze(sum(bsxfun(@times,pair_vm,Qbar(:,:,pair)),1)).*squeeze(pair_vm),1)';
    
    fprintf('%d\n', i);
end

clf
drawMesh(V, F, 'facecolor','y', 'edgecolor','k', 'linewidth', 1.2);
view([0 90])
axis equal
axis off
camlight
lighting gouraud
cameramenu