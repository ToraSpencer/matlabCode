


delter_near_slice_centre_point = 0;
delter_near_slice_centre_point1 = cell(1);
P_slice_zero_one_matrix = cellfun('isempty',P_slice_centre);
[M_P_slice_center,N_P_slice_center] = size(P_slice_zero_one_matrix);

for i = 1:num_img
    i
    if(i == 1)
        for j = 1:N_P_slice_center-sum(P_slice_zero_one_matrix(i,:),2)
            P_uv_data_related{j,:} = P_u_sampledot{i,j};
        end
    else
        for n = 1:N_P_slice_center-sum(P_slice_zero_one_matrix(i,:),2)
            for m = 1:N_P_slice_center-sum(P_slice_zero_one_matrix(i-1,:),2)
                delter_near_slice_centre_point1{i,n}(:,m) = sqrt(...
                    (P_slice_centre{i,n}(:,1) - P_slice_centre{i-1,m}(:,1))^2 ...
                    + (P_slice_centre{i,n}(:,2) - P_slice_centre{i-1,m}(:,2))^2);
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
P_uv_data_related_trans = P_uv_data_related.';
P_uv_data_related_zero_one_matrix_trans_trans = cellfun('isempty',P_uv_data_related_trans);
[M_final_trans_P_uv_data_related,N_trans_final_P_uv_data_related] = size(P_uv_data_related_zero_one_matrix_trans_trans);

N_trans_final_P_uv_data_related = 2;

for k = 1:N_trans_final_P_uv_data_related
    N_effective = M_final_trans_P_uv_data_related - sum(P_uv_data_related_zero_one_matrix_trans_trans(:,k));
    if(N_effective <= 3)
        continue;
    end
    for j = 1:10
        for i = 1:N_effective
            P_uv_data{k,j}(:,i) = P_uv_data_related_trans{i,k}(:,j);
        end
            count_v = count_v +1;
            trans_P_uv_data = P_uv_data{k,j}.';
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
            sp_v = spmak(NodeVector_v,trans_P_v);
            fnplt(sp_v,[NodeVector_v(1),NodeVector_v(end)]);
    end
end

piece_u = 100;   %就有30个横坐标上的值
piece_v = 100;
u = linspace(0, 1-0.0001, piece_u);
X_u = zeros(1, piece_u);
Y_u = zeros(1, piece_u);
Z_u = zeros(1, piece_u);
P_surface_u_slice = cell(1);
for k = 1:N_trans_final_P_uv_data_related
    for j = 1:M_final_trans_P_uv_data_related - sum(P_uv_data_related_zero_one_matrix_trans_trans(:,k))
        for i = 1:10
            %         P_surface_u_slice得到每一行坐标
            %         x
            %         y
            %         z
            P_surface_u_slice{k,j}(:,i) = P_surface{k,i}(:,j);
        end
    end
end
for k = 1:N_trans_final_P_uv_data_related
    for i = 1:M_final_trans_P_uv_data_related - sum(P_uv_data_related_zero_one_matrix_trans_trans(:,k))   %i：第几层轮廓
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
            for ii = 0 : 1 : M_final_trans_P_uv_data_related - sum(P_uv_data_related_zero_one_matrix_trans_trans(:,k))-1
                Nik_v(ii+1, 1) = BaseFunction(ii, 3, v(j), v_sampledot_vector{2});
            end
            X_uv(j, i) = Nik_v' * X_u(:, i);
            Y_uv(j, i) = Nik_v' * Y_u(:, i);
            Z_uv(j, i) = Nik_v' * Z_u(:, i);
        end
    end
    mesh(X_uv,Y_uv,Z_uv);
end
%}
xlabel('X');
ylabel('Y');
zlabel('Z');
view(90, -45);
toc;