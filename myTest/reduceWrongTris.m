function [tris ] = reduceWrongTris( inputTris )

% 去除面信息中错误的三角片（同一个三角片中存在两个相同的点索引）

index1 = inputTris(:, 1); 
index2 = inputTris(:, 2); 
index3 = inputTris(:, 3); 
compare12 = (index1 == index2);
compare23 = (index2 == index3);
compare13 = (index1 == index3);
compare = (compare12 | compare13 | compare23);
compare = ~compare;
tris = inputTris(compare, :);

end

