function sampledot = SampleFromOutline(Q,num_single_img_outline)
sampledot = zeros(1,num_single_img_outline);
num_data_dot = length(Q);

if((4 <= num_data_dot) && (num_data_dot <= 8))
    sampledot = Q;
end
if(9 <= num_data_dot) && (num_data_dot <= 11)
    length_bi = floor(num_data_dot/2);
    for k = 1:length_bi
        sampledot(k,:) = Q(2*k,:);
    end
end
if(12 <= num_data_dot) && (num_data_dot <= 20)
    length_bi = floor(num_data_dot/3);
    for k = 1:length_bi
        sampledot(k,:) = Q(3*k,:);
    end
end
if(num_data_dot > 20)
    length_bi = floor(num_data_dot/5);
    for k = 1:length_bi
        sampledot(k,:) = Q(5*k,:);
    end
end
