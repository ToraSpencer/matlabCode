function [in] = segmentation_region_grow(v,f,points)
%该程序用于在知道边界点的情况下，将网格沿着边界点分为上下两部分
% v - 需要被分割的体的点坐标
% f - 需要被分割的体的面拓扑
% points - 边界线的点坐标
% by lsp
%找最近的两个点，并给出其对应的序号,将与其相邻的序号考虑删除
% center = mean(v);
for i = 1:length(points)
    
    d = sqrt(sum((points(i,:)-v).^2,2));
    ii = find(d == min(d));
%     d(ii) = max(d);
%     jj = find(d == min(d));
%     In(i,:) = [ii,jj];

    IN(i) = ii(1);
    
end
edges = [sort([ f(:,1)  f(:,2) ]')';sort([ f(:,1)  f(:,3) ]')';sort([ f(:,2)  f(:,3) ]')'];
%去掉重边
[etemp1,~,~]=unique(edges,'rows','stable');
%找到牙龈线在的边
[ia2,ib]=find(bsxfun(@eq,etemp1(:,1),IN));
a = etemp1(ia2,:);
[ia3,ib3]=find(bsxfun(@eq,a(:,2),IN));
b = a(ia3,:);%b为牙龈线所在的边的集合
% [ia2,ib2]=find(bsxfun(@eq,etemp1(:,2),IN));
% a2 = etemp1(ia2,:);
% A = [a;a2];
%%用区域增长法，将两边分隔开
%1.随机选一个面片，以该面片为中心（选取第一个面片）
ff = f;neblist = [];
seed = ff(1,:);
%进入循环
E = [sort([seed(:,1) seed(:,2)]')';sort([seed(:,2) seed(:,3)]')';sort([seed(:,1) seed(:,3)]')'];
neblist = [neblist;seed];


ff(1,:) = 0;

outline  = b;
KK = [];
% figure()
% trimesh(f,v(:,1),v(:,2),v(:,3))
% hold on
while  any(outline(:,1)) && any(ff(:,1)) && ~isempty(E)
    [c, iaed, ibed] = intersect(b,E,'rows');
%     plot_edges(v,E,'r','LineWidth',2);
    if ~isempty(c)

        %检查是否遍历完，如未遍历结束循环继续,并且与牙龈线边界贴合的边，不做考虑
        outline(iaed,1) = 0;
        E(ibed,:) = [];
    end
    l = length(E);
    E1 = sort([ ff(:,1)  ff(:,2) ]')';
    E2 = sort([ ff(:,2)  ff(:,3) ]')';
    E3 = sort([ ff(:,1)  ff(:,3) ]')';
    K = [];
    
%     for i = 1:l
%         H1 = ismember(E1,E(i,:),'rows');
%         H2 = ismember(E2,E(i,:),'rows');
%         H3 = ismember(E3,E(i,:),'rows');
%         k1 = find( H1 == 1 );k2 = find( H2 == 1 );k3 = find( H3 == 1 );
%         K = [K;k1;k2;k3];    
%     end
    [ed1, k1, ~] = intersect(E1,E,'rows');[ed2, k2, ~] = intersect(E2,E,'rows');[ed3, k3, ~] = intersect(E3,E,'rows');
%     E = [ed1;ed2;ed3];
    K = [K;k1;k2;k3];  
    %去掉重复的K
    K = unique(K);
    
    %去掉与原先相重复的
%     [c1, ia1, ib1] = intersect(KK,K');
%     if isempty(c1)
%        KK = [KK;K];
%     else
%        K(ib1) = [];
%        KK = [KK;K];
%     end 
    k = find(ff(K,1) ~= 0);
    K = K(k);
    KK = [KK;K];
    seed = f(K,:);
    ff(K,1) = 0;
%     E = [[seed(:,1) seed(:,2)];[seed(:,2) seed(:,3)];[seed(:,1) seed(:,3)]];
    E = [sort([seed(:,1) seed(:,2)]')';sort([seed(:,2) seed(:,3)]')';sort([seed(:,1) seed(:,3)]')'];
    %判断增长出来的边是否触碰到了牙龈线边界
%     [c, iaed, ibed] = intersect(b,E,'rows');
%     plot_edges(v,E,'r','LineWidth',2);
    
%     if ~isempty(c)
% 
%         %检查是否遍历完，如未遍历结束循环继续,并且与牙龈线边界贴合的边，不做考虑
%         outline(iaed,1) = 0;
%         E(ibed,:) = [];
%     end
    
end
KK = unique(KK);
    
    







 in = KK;












