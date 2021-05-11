clc;
clear;
tic;
tifFile=dir(['E:/tif_mask_original' ,'/*.tif']);

num_img = 5;%length(tifFile);
is_uvector_count = true;
sample_spline_points = zeros(3,10);
is_vvector_count = true;
slice_centre_dot = zeros(25,3);
next_slice_centre_dot = zeros(25,3);
P_u_sampledot = cell(1);
P_slice_centre = cell(1);
u_sampledot_vector = cell(1,2);
P_uv_data_related = cell(1);
P_uv_data = cell(1);
P_surface = cell(1,10);
v_sampledot_vector = cell(1,2);
for i = 1:num_img
    start_offset = 33;
    file_name=['E:/tif_mask_original' '/' tifFile(i).name];
    I = imread(file_name);
    bw = im2bw(I);
    [B,L] = bwboundaries(bw);
    del_foreground = I*0;
    axis on;
    num_single_img_outline = length(B);
    %     fprintf('The  117608_0%3.3d num_single_img_outline is %8.3d\n',start_offset+i,num_single_img_outline);
    sampledot = cell(1,num_single_img_outline);
    for j = 1: num_single_img_outline
        
        num_data_dot = length(B{j});
        if(num_data_dot <= 3)
            B{j} = 0;
            continue;
        end
        if((4 <= num_data_dot) && (num_data_dot <= 8))
            sampledot{j} = B{j};
        end
        if(9 <= num_data_dot) && (num_data_dot <= 11)
            length_bi = floor(num_data_dot/2);
            for k = 1:length_bi
                sampledot{j}(k,:) = B{j}(2*k,:);
            end
        end
        if(12 <= num_data_dot) && (num_data_dot <= 20)
            length_bi = floor(num_data_dot/3);
            for k = 1:length_bi
                sampledot{j}(k,:) = B{j}(3*k,:);
            end
        end
        if(num_data_dot > 20)
            length_bi = floor(num_data_dot/5);
            for k = 1:length_bi
                sampledot{j}(k,:) = B{j}(5*k,:);
            end
        end
        length_sampledot = length(sampledot{j});
        if(~((sampledot{j}(1,1) == sampledot{j}(length_sampledot,1,1)) ...
                && (sampledot{j}(1,2) == sampledot{j}(length_sampledot,2))))
            sampledot{j}(length_sampledot+1,:) = sampledot{j}(1,:);
        end
        l_sampledot = length(sampledot{j});
        for k = 1:l_sampledot
            sampledot{j}(k,3) = i + start_offset;
        end
        hold on;
        PreNodeVector_u = PreChordLength_Para(sampledot{j});
        NodeVector_u = ChordLength_Para(sampledot{j},PreNodeVector_u,3);
        P = GetControlMatrix(sampledot{j},PreNodeVector_u,NodeVector_u);
        sp = spmak(NodeVector_u,P.');
        %         fnplt(sp,[NodeVector_u(1),NodeVector_u(end)]);
        spline_points = fnplt(sp,[NodeVector_u(1),NodeVector_u(end)]);
        num_spline_points = length(spline_points);
        inter = floor((num_spline_points)/9);
        %在样条曲线上取9个采样点,其中一个点是重合的
        for k = 1:9
            sample_spline_points(:,k) = spline_points(:,inter*k);
        end
        sample_spline_points(:,10) = sample_spline_points(:,1);
        %得出10个控制点
        if(is_uvector_count == true)
            u_surface_prenode =  PreChordLength_Para(sample_spline_points.');
            u_surface_node = ChordLength_Para(sample_spline_points.',u_surface_prenode,3);
            u_P = GetControlMatrix(sample_spline_points.',u_surface_prenode,u_surface_node);
            u_sampledot_vector{1} = u_surface_prenode;
            u_sampledot_vector{2}=  u_surface_node;
            is_uvector_count = false;
        end
        if(is_uvector_count == false)
            u_P = GetControlMatrix(sample_spline_points.',u_sampledot_vector{1},u_sampledot_vector{2});
        end
        X_mean = sum(sample_spline_points(1,:)')/10;
        Y_mean = sum(sample_spline_points(2,:)')/10;
        Z_mean = sum(sample_spline_points(3,:)')/10;
        P_slice_centre{i,j} = [X_mean Y_mean Z_mean];
        %         for t = 0:pi/3:(2*pi)
        %             C = Y_mean - tand(t)*X_mean;
        %             for p = 1:num_spline_points
        %                 d = abs((tand(t)*spline_points(1,p) - spline_points(2,p) + C)/(tand(t)^2+1));
        % %                 line(C,spline_points());
        %             end
        %
        %         end
        trans_u_P = u_P.';
        %保留第i层第j个轮廓的10个u向的控制点
        P_u_sampledot{i,j} = trans_u_P;
        %         plot3(trans_u_P(1,:),trans_u_P(2,:),trans_u_P(3,:),'r+');
        %         plot3(sample_spline_points(1,:),sample_spline_points(2,:),sample_spline_points(3,:),'r*');
        %{
        if(mod(i,2) == 1)
            slice_centre_dot(j,:) = sum(sample_spline_points.')/10;
            %           把行向量为0的行去除掉
            slice_centre_dot(all(slice_centre_dot == 0,2),:) = [];
            
            %             slice_centre_dot(j,3) = slice_centre_dot(j,3) + 1
        end
        if(mod(i,2) == 0)
            next_slice_centre_dot(j,:) = sum(sample_spline_points.')/10;
            next_slice_centre_dot(all(next_slice_centre_dot == 0,2),:) = [];
            
        end
        %}
        
        %         for t=0:pi/10:2*pi
        %             m = (Y_mean - tand(t)*X_mean)*(sum(P(:,2))-sum(P(:,1)))
        %             sample_spline_points1 = m*P(:,1:2)
        %             plot(sample_spline_points1(:,1),sample_spline_points1(:,2),'ro');
        %         end
        %         fprintf('The  X_mean is %3.3d Y_mea nis is %3.3d Z_mean is %3.3d\n',X_mean,Y_mean,Z_mean);
        %         plot3(X_mean,Y_mean,Z_mean,'r*');
        %         plot(X_mean,Y_mean,'r*');
        %                 fill(points(1,:),points(2,:),'w');
       
        %         P_uv_data =
        %         2×1 cell 数组
        %         {3×10 double}
        %         {3×10 double}
        %         if(i == 1)
        %             P_uv_data{j,i} = P_u_sampledot{i,j};
        %         end
    end
     %{
    if(i > 1)
        i
        num_slice_centre_dot = size(slice_centre_dot,1)
        num_next_slice_centre_dot = size(next_slice_centre_dot,1)
        delter_near_slice_centre_point = zeros(num_slice_centre_dot,num_next_slice_centre_dot);
        for m = 1:num_slice_centre_dot
            for n = 1:num_next_slice_centre_dot
                delter_near_slice_centre_point(n,m) = sqrt(...
                    (slice_centre_dot(m,1) - next_slice_centre_dot(n,1))^2 ...
                    + (slice_centre_dot(m,2) - next_slice_centre_dot(n,2))^2);
                if(delter_near_slice_centre_point < 30)
%                     P_uv_data(:,i) = P_u_sampledot{i,1}(:,j);
                end
            end
        end
        if(mod(i,2) == 1)
            next_slice_centre_dot = [];
            num_next_slice_centre_dot = 0;
        end
        if(mod(i,2) == 0)
            slice_centre_dot = [];
            num_slice_centre_dot = 0;
        end
    end
        %}
        hold on;
end


delter_near_slice_centre_point = 0;
% num_nonzero_cell = 0;
P_slice_zero_one_matrix = cellfun('isempty',P_slice_centre);
[M_P_slice_center,N_P_slice_center] = size(P_slice_zero_one_matrix);
% nonzero_P_centre_matrix = (P_slice_zero_one_matrix==0);
for i = 1:num_img
    i
    if(i == 1)
        for j = 1:N_P_slice_center-sum(P_slice_zero_one_matrix(i,:),2)
            P_uv_data_related{j,:} = P_u_sampledot{i,j};
        end
    else
        for n = 1:N_P_slice_center-sum(P_slice_zero_one_matrix(i,:),2)
            for m = 1:N_P_slice_center-sum(P_slice_zero_one_matrix(i-1,:),2)
                delter_near_slice_centre_point = sqrt(...
                    (P_slice_centre{i,n}(:,1) - P_slice_centre{i-1,m}(:,1))^2 ...
                    + (P_slice_centre{i,n}(:,2) - P_slice_centre{i-1,m}(:,2))^2);
                if(delter_near_slice_centre_point <= 5)
                    issave = true;
                    P_uv_data_zero_one_metrix = cellfun('isempty',P_uv_data_related);
                    nonzero_m_zero_one_matrix = (P_uv_data_zero_one_metrix == 0);
                    count_row_nonzero_P_uv_data = sum(nonzero_m_zero_one_matrix(m,:),2);
                    P_uv_data_related{m,(count_row_nonzero_P_uv_data+1)} = P_u_sampledot{i,n};
                    break;
                end
            end
            if(delter_near_slice_centre_point > 5)
               P_uv_data_zero_one_metrix = cellfun('isempty',P_uv_data_related);
               [M_P_uv_data_related,N_P_uv_data_related] = size(P_uv_data_zero_one_metrix);
               P_uv_data_related{M_P_uv_data_related+1,1} = P_u_sampledot{i,n};
            end
        end
    end
end

count_v = 0;
trans_P_uv_data_related = P_uv_data_related.';
trans_P_uv_data_related_zero_one_matrix = cellfun('isempty',trans_P_uv_data_related);
[M_final_trans_P_uv_data_related,N_trans_final_P_uv_data_related] = size(trans_P_uv_data_related_zero_one_matrix);

for k = 1:N_trans_final_P_uv_data_related
    for j = 1:10
        for i = 1:M_final_trans_P_uv_data_related - sum(trans_P_uv_data_related_zero_one_matrix(:,k))
            P_uv_data{k,j}(:,i) = trans_P_uv_data_related{i,k}(:,j);
        end
        if(k <= 2)
            count_v = count_v +1;
            trans_P_uv_data = P_uv_data{k,j}.'
            if(is_vvector_count == true)
                PreNodeVector_v = PreChordLength_Para(trans_P_uv_data);
                NodeVector_v = ChordLength_Para(trans_P_uv_data,PreNodeVector_v,3);
                P_v = GetControlMatrix(trans_P_uv_data,PreNodeVector_v,NodeVector_v);
                v_sampledot_vector{1} = PreNodeVector_v;
                v_sampledot_vector{2} = NodeVector_v;
                is_vvector_count = false;
            end
            if(is_vvector_count == false)
                P_v = GetControlMatrix(trans_P_uv_data, v_sampledot_vector{1},v_sampledot_vector{2});
            end
            trans_P_v = P_v.';
            %在u方向的控制点取纵向的点，得出整个曲面的控制点
            P_surface{k,j} = trans_P_v;
            trans_P_uv_data = 0;
            %     plot3(trans_P_v(1,:),trans_P_v(2,:),trans_P_v(3,:),'square');
            %     Length_v = length(NodeVector_v);
            %     sp_v = spmak(NodeVector_v,trans_P_v);
            %     fnplt(sp_v,[NodeVector_v(1),NodeVector_v(Length_v)]);
        end
    end
end

%{
for j = 1:10
    for i = 1:num_img
%把
        P_uv_data{j}(:,i) = P_u_sampledot{i,1}(:,j);
    end
    count_v = count_v +1;
    trans_P_uv_data = P_uv_data{j}.';
    if(is_vvector_count == true)
        PreNodeVector_v = PreChordLength_Para(trans_P_uv_data);
        NodeVector_v = ChordLength_Para(trans_P_uv_data,PreNodeVector_v,3);
        P_v = GetControlMatrix(trans_P_uv_data,PreNodeVector_v,NodeVector_v);
        v_sampledot_vector{1} = PreNodeVector_v;
        v_sampledot_vector{2} = NodeVector_v;
        is_vvector_count = false;
    end
    if(is_vvector_count == false)
        P_v = GetControlMatrix(trans_P_uv_data, v_sampledot_vector{1},v_sampledot_vector{2});
    end
    trans_P_v = P_v.';
    %在u方向的控制点取纵向的点，得出整个曲面的控制点
    P_surface{j} = trans_P_v;
%     plot3(trans_P_v(1,:),trans_P_v(2,:),trans_P_v(3,:),'square');
%     Length_v = length(NodeVector_v);
%     sp_v = spmak(NodeVector_v,trans_P_v);
    %     fnplt(sp_v,[NodeVector_v(1),NodeVector_v(Length_v)]);
end
%}
piece_u = 30;   %就有30个横坐标上的值
piece_v = 30;
u = linspace(0, 1-0.0001, piece_u);
X_u = zeros(1, piece_u);
Y_u = zeros(1, piece_u);
Z_u = zeros(1, piece_u);
P_surface_u_slice = cell(1,5);
for k = 1:2
    for j = 1:5
        for i = 1:10
            %         P_surface_u_slice得到每一行坐标
            %         x
            %         y
            %         z
            P_surface_u_slice{k,j}(:,i) = P_surface{k,i}(:,j);
        end
    end
end
for k = 1:2
    for i = 1:5   %i：第几层轮廓
        for j = 1 : piece_u
            for ii = 0 : 1: 9
                Nik_u(ii+1, 1) = BaseFunction(ii, 3 , u(j), u_sampledot_vector{2});
            end
            X_u(i,j) = P_surface_u_slice{k,i}(1,:)*Nik_u;
            Y_u(i,j) = P_surface_u_slice{k,i}(2,:)*Nik_u;
            Z_u(i,j) = P_surface_u_slice{k,i}(3,:)*Nik_u;
        end
    end
    X_uv = zeros(piece_v, piece_u);
    Y_uv = zeros(piece_v, piece_u);
    Z_uv = zeros(piece_v, piece_u);
    v = linspace(0, 1-0.0001, piece_v);
    for i = 1 : piece_u
        for j = 1 : piece_v
            for ii = 0 : 1 : 4
                Nik_v(ii+1, 1) = BaseFunction(ii, 3, v(j), v_sampledot_vector{2});
            end
            X_uv(j, i) = Nik_v' * X_u(:, i);
            Y_uv(j, i) = Nik_v' * Y_u(:, i);
            Z_uv(j, i) = Nik_v' * Z_u(:, i);
        end
    end
    mesh(X_uv,Y_uv,Z_uv);
    % plot3(X_uv,Y_uv,Z_uv)
end
%}
xlabel('X');
ylabel('Y');
zlabel('Z');
view(90, -45);
toc;