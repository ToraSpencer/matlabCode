clc;
clear all;
tic;

tifFile = dir(['E:/tif_mask_original_after_segment' ,'/*.tif']);
num_img = length(tifFile);

initial_image = 92;
count = 0;
for i = 1:num_img
    if(mod(i,2) == 0)
        count = count + 1;
        start_offset = 33;
        current_i = i - start_offset;
        filename = ['E:/tif_mask_original_after_segment/' tifFile(i).name];
        
        G = imread(filename);
        I = im2bw(double(G));
%         bwmorph3
        bw_I = imopen(I,ones(9,9));
        F = getframe(gca);
        
%         if ~isdir('E:/test_Bspline_layer_tif')
%             mkdir('E:/test_Bspline_layer_tif');
%         end
%         
        if ~isdir('E:/test_segment_layer_2_bw_I_tif')
            mkdir('E:/test_segment_layer_2_bw_I_tif');
        end
        
        imwrite(bw_I,['E:/test_segment_layer_2_bw_I_tif/117608_0',num2str(count,'%3.3d'),'.tif'],'tif');
%         imwrite(I,['E:/test_Bspline_layer_tif/117608_0',num2str(count,'%3.3d'),'.tif'],'tif');
    end
end