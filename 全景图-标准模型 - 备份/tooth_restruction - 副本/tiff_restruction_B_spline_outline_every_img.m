clc;
clear all;

% tifFile=dir(['E:/tif_test1' ,'/*.tif']);
tifFile=dir(['E:/tif_mask_original' ,'/*.tif']);
num_img = length(tifFile);

AllFile = cell(num_img,2);
rotate_matrix = [0 1;-1 0];
count = 0;

if ~isdir('E:/tif_segment_Bspline')
    mkdir('E:/tif_segment_Bspline');
end

for i = 1:10
    start_offset = 33;
    AllFile(i,1) = {i+start_offset};
    file_name=['E:/tif_mask_original' '/' tifFile(i).name];
    I = imread(file_name);
    bw = im2bw(I);
    img_filled_holes = imfill(bw,'holes');
    %     imshow(img_filled_holes);
    [B,L] = bwboundaries(img_filled_holes);
    del_foreground = I*0;
    img_orig = imshow(del_foreground,'border','tight');
    hold on;
    num_single_img_outline = length(B);
%     fprintf('The  117608_0%3.3d num_single_img_outline is %8.3d\n',start_offset+i,num_single_img_outline);
    sampledot = cell(1,num_single_img_outline);
%         figure('Color','black');
    %     title(['order=',num2str(start_offset+i)])
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
        if(~((sampledot{j}(1,1) == sampledot{j}(length_sampledot,1)) ...
                && (sampledot{j}(1,2) == sampledot{j}(length_sampledot,2))))
            sampledot{j}(length_sampledot+1,:) = sampledot{j}(1,:);
        end
        %         hold on;
        axis off;
        %         plot3(sampledot{j}(:,1),sampledot{j}(:,2),sampledot{j}(:,3),'ro');
        %         trans_sampledot = rotate_matrix*sampledot{j}.';
        %         plot(trans_sampledot(1,:),trans_sampledot(2,:),'ro');
        PreNodeVector_u = PreChordLength_Para(sampledot{j});
        NodeVector_u = ChordLength_Para(sampledot{j},PreNodeVector_u,3);
        D = GetCoefMatrix(sampledot{j},PreNodeVector_u,NodeVector_u);
        display_matrix = rotate_matrix*(D.');
        %         plot(D(:,1),D(:,2),'r+');
        %         plot(D(:,1),D(:,2));
        Length_u = length(NodeVector_u);
        sp = spmak(NodeVector_u,D.');
        fnplt(sp,[NodeVector_u(1),NodeVector_u(Length_u)]);
        points = fnplt(sp,[NodeVector_u(1),NodeVector_u(Length_u)]);
        count = count + 1;
        fprintf('The points is %8d count is %4d\n',length(points),count);
        fill(points(1,:),points(2,:),'w');
    end

    data = getframe(gcf);
    gffim = frame2im(data);
%     imwrite(data.cdata.*255,['E:/tif_test2/117608_0',num2str(start_offset+i,'%3.3d'),'.tif'],'tif');
    imwrite(data.cdata.*255,['E:/tif_segment_Bspline/117608_0',num2str(start_offset+i,'%3.3d'),'.tif'],'tif');
    pause(0.3);
end
%
xlabel('X');
ylabel('Y');
zlabel('Z');
view(90, -45);