clc;
clear all;
tic;
tifFile=dir(['E:/tif_mask_original' ,'/*.tif']);
num_img = length(tifFile);

AllFile = cell(num_img,2);
count = 0;
slice_centre_dot = zeros(25,3);
sample_spline_points = zeros(3,10);
next_slice_centre_dot = zeros(25,3);
next_sample_spline_points = zeros(3,10);
P_u_sampledot = cell(4,2);
u_sampledot_vector = cell(1,2);
for i = 1:5
    count = count +1
    start_offset = 33;
    AllFile(i,1) = {i+start_offset};
    file_name=['E:/tif_mask_original' '/' tifFile(i).name];
    I = imread(file_name);
    bw = im2bw(I);
    %     img_filled_holes = imfill(bw,'holes');
    %         imshow(I);
    [B,L] = bwboundaries(bw);
    del_foreground = I*0;
    %     img_orig = imshow(del_foreground,'border','tight');
    %     figure;
    %     hold on;
    axis on;
    num_single_img_outline = length(B);
    fprintf('The  117608_0%3.3d num_single_img_outline is %8.3d\n',start_offset+i,num_single_img_outline);
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
        %         figure;
        hold on;
        %         plot(sampledot{j}(:,2),sampledot{j}(:,1),'ro');
        PreNodeVector_u = PreChordLength_Para(sampledot{j});
        NodeVector_u = ChordLength_Para(sampledot{j},PreNodeVector_u,3);
        P = GetControlMatrix(sampledot{j},PreNodeVector_u,NodeVector_u);
        %         display_matrix = rotate_matrix*(D.');
        Length_u = length(NodeVector_u);
        sp = spmak(NodeVector_u,P.');
%         fnplt(sp,[NodeVector_u(1),NodeVector_u(Length_u)])
        spline_points = fnplt(sp,[NodeVector_u(1),NodeVector_u(Length_u)]);
        num_spline_points = length(spline_points);
        inter = floor((num_spline_points)/9);
        %在样条曲线上取9个采样点,其中一个点是重合的
        for k = 1:9
            sample_spline_points(:,k) = spline_points(:,inter*k);
        end
        sample_spline_points(:,10) = sample_spline_points(:,1);
        %得出10个控制点
        if(count <= 1)
            u_surface_prenode =  PreChordLength_Para(sample_spline_points.');
            u_surface_node = ChordLength_Para(sample_spline_points.',u_surface_prenode,3);
            u_P = GetControlMatrix(sample_spline_points.',u_surface_prenode,u_surface_node);
            u_sampledot_vector{1} = u_surface_prenode;
            u_sampledot_vector{2}=  u_surface_node;
        end
        if(count >1)
            u_P = GetControlMatrix(sample_spline_points.',u_sampledot_vector{1},u_sampledot_vector{2});
        end

        trans_u_P = u_P.';
        %保留每一层的控制点
        P_u_sampledot{i,j} = trans_u_P;
%         plot3(trans_u_P(1,:),trans_u_P(2,:),trans_u_P(3,:),'r+');
%         plot3(sample_spline_points(1,:),sample_spline_points(2,:),sample_spline_points(3,:),'r*');
        if(mod(i,2) == 1)
            slice_centre_dot(j,:) = sum(sample_spline_points.')/10;
            slice_centre_dot(all(slice_centre_dot == 0,2),:) = [];
            slice_centre_dot(j,3) = slice_centre_dot(j,3) + 1;
        end
        if(mod(i,2) == 0)
            next_slice_centre_dot(j,:) = sum(sample_spline_points.')/10;
            next_slice_centre_dot(all(next_slice_centre_dot == 0,2),:) = [];
        end
        X_mean = sum(sample_spline_points(1,:)')/10;
        Y_mean = sum(sample_spline_points(2,:)')/10;
        Z_mean = sum(sample_spline_points(3,:)')/10;
        %         for t=0:pi/10:2*pi
        %             m = (Y_mean - tand(t)*X_mean)*(sum(P(:,2))-sum(P(:,1)))
        %             sample_spline_points1 = m*P(:,1:2)
        %             plot(sample_spline_points1(:,1),sample_spline_points1(:,2),'ro');
        %         end
        %         fprintf('The  X_mean is %3.3d Y_meanis is %3.3d Z_mean is %3.3d\n',X_mean,Y_mean,Z_mean);
        %         plot3(X_mean,Y_mean,Z_mean,'r*');
        %         plot(X_mean,Y_mean,'r*');
        %                 fill(points(1,:),points(2,:),'w');
    end
 
    if(i > 1)
        num_slice_centre_dot = size(slice_centre_dot,1);
        num_next_slice_centre_dot = size(next_slice_centre_dot,1);
        delter_near_slice_centre_point = zeros(num_slice_centre_dot,num_next_slice_centre_dot);
        for m = 1:num_slice_centre_dot
            for n = 1:num_next_slice_centre_dot
                delter_near_slice_centre_point(n,m) = sqrt(...
                    (slice_centre_dot(m,1) - next_slice_centre_dot(n,1))^2 ...
                    + (slice_centre_dot(m,2) - next_slice_centre_dot(n,2))^2);
            end
            if(delter_near_slice_centre_point < 10)
                
            end
        end
    end
    
    hold on;
end

P_uv_data = cell(1,10);
P_surface = cell(1,10);
v_sampledot_vector = cell(1,2);
count_v = 0;
for j = 1:10
    for i = 1:5
        P_uv_data{j}(:,i) = P_u_sampledot{i,1}(:,j);
    end
    count_v = count_v +1;
    trans_P_uv_data = P_uv_data{j}.';
    if(count_v <=1)
        PreNodeVector_v = PreChordLength_Para(trans_P_uv_data);
        NodeVector_v = ChordLength_Para(trans_P_uv_data,PreNodeVector_v,3);
        P_v = GetControlMatrix(trans_P_uv_data,PreNodeVector_v,NodeVector_v);
        v_sampledot_vector{1} = PreNodeVector_v;
        v_sampledot_vector{2} = NodeVector_v;
        trans_P_v = P_v.';
    end
    if(count_v > 1)
        P_v = GetControlMatrix(trans_P_uv_data, v_sampledot_vector{1},v_sampledot_vector{2});
    end
    trans_P_v = P_v.';
    %在u方向的控制点每圈取5个，得出整个曲面的控制点
    P_surface{j} = trans_P_v;
%     plot3(trans_P_v(1,:),trans_P_v(2,:),trans_P_v(3,:));
%     plot3(trans_P_v(1,:),trans_P_v(2,:),trans_P_v(3,:),'square');
    Length_v = length(NodeVector_v);
    sp_v = spmak(NodeVector_v,trans_P_v);
%     fnplt(sp_v,[NodeVector_v(1),NodeVector_v(Length_v)]);
end

piece_u = 30;   %就有30个横坐标上的值
piece_v = 30;
u = linspace(0, 1-0.00001, piece_u);
X_u = zeros(1, piece_u);
Y_u = zeros(1, piece_u);
Z_u = zeros(1, piece_u);
P_surface_u_slice = cell(1,5);
for j = 1:5
    for i = 1:10
        P_surface_u_slice{j}(:,i) = P_surface{i}(:,j);
    end
end

for i = 1:5   %i：第几层轮廓
    for j = 1 : piece_u
        for ii = 0 : 1: 9
            Nik_u(ii+1, 1) = BaseFunction(ii, 3 , u(j), u_sampledot_vector{2})
            fprintf('The Nik_u is %8.5d\n',Nik_u);
        end
        X_u(i,j) = P_surface_u_slice{i}(1,:)*Nik_u;
        Y_u(i,j) = P_surface_u_slice{i}(2,:)*Nik_u;
        Z_u(i,j) = P_surface_u_slice{i}(3,:)*Nik_u;
    end
end

X_uv = zeros(piece_v, piece_u);
Y_uv = zeros(piece_v, piece_u);
Z_uv = zeros(piece_v, piece_u);
v = linspace(0.01, 1-0.0001, piece_v);
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

xlabel('X');
ylabel('Y');
zlabel('Z');
view(90, -45);
%}
toc;