clc;
clear all;
close all;

%% 读取数据
load('data.mat');

        [Vj, Fj] = readOBJ('lower_jaw.obj');
 
        jawPts = [
            -32.0, -36.0, 0
            -28.5, -27.5, 0
            -24.5, -20.0, 0
            -21.0, -13.0, 0
            -15.5, -7.0, 0
            -10.0, -2.0, 0
            -4.0, -0.5, 0
            4.0, -0.5, 0
            10.0, -2.0, 0
            15.5, -7.0, 0
            21.0, -13.0, 0
            24.5, -20.0, 0
            28.5, -27.5, 0
            32.0, -36.0, 0];
        jawPts = jawPts(end:-1:1,:);
 
    
    dentalCenter = handle.dentalCenter;
    dentalFrame = handle.dentalFrame;
	toothdata = handle.model.toothdata;
    nT = size(toothdata,2);
 
    % 牙齿中心
    toothCenter = zeros(nT,3);
    for i = 1:nT
        toothCenter(i,:) = mean(toothdata{1,i}(toothdata{3,i},:));
    end
    
    
%%
    % 将牙齿中心转换到牙颌xy平面上并用标准椭圆拟合
    pts0 = bsxfun(@minus, toothCenter, dentalCenter) * dentalFrame;
    coff = fit_ellipse(pts0(:,1), pts0(:,2), 'standard');
    ax_o = [-coff(3)/(2*coff(1)), -coff(4)/(2*coff(2))];
    temp = coff(3)*coff(3)/(4*coff(1)) + coff(4)*coff(4)/(4*coff(2)) - coff(5);
    ax_l = sqrt(temp/coff(1));
    ax_s = sqrt(temp/coff(2));
    
    % 均匀选取两端牙齿之间的点
    step = (-pi/2:0.01:3*pi/2)';
    pts = [ax_l*cos(step), ax_s*sin(step)];
    pts = bsxfun(@plus, pts, ax_o);
    idx1 = knnsearch(pts, (pts0(1,1:2) + [-pts0(nT,1), pts0(nT,2)])/2);
    idx2 = knnsearch(step, pi-step(idx1));
    sel = round(linspace(idx1, idx2, 14));
    ptsn = pts(sel,:);
  
    % 计算权重
    weight = bsxfun(@minus, Vj(:,1:2), reshape(jawPts(:,1:2)', [1 2 14]));
    weight = sum(weight.^2, 2)/1000;
    weight = exp(-weight);
    weight = bsxfun(@rdivide, weight, sum(weight,3));
  
    % 标准颌骨变形
    Vj(:,1:2) = sum( bsxfun(@times, reshape((ptsn - jawPts(:,1:2))', 1, 2, nT), weight), 3 ) + Vj(:,1:2);
    
    % 转换到全局坐标系
	Vj = bsxfun(@plus, Vj * dentalFrame', dentalCenter);
    
	Vall = [];
    Fall = [];
    for i = 1:nT
        Fall = [Fall; toothdata{10,i} + size(Vall,1)];
        Vall = [Vall; toothdata{9,i}];
    end
    
    
    %%
    figure
    drawMesh(Vall, Fall, 'facecolor','y', 'edgecolor','none', 'facealpha',1.0);
    drawMesh(Vj, Fj, 'facecolor','g', 'edgecolor','none', 'facealpha',0.8);

    view(3)
    axis equal
    axis off
    camlight
    lighting gouraud
    set(gca, 'Position',[0 0 1 1]);
 
 