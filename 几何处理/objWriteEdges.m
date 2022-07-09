function objWriteEdges( fileName, edges, meshVers)
% 输出表示边的OBJ文件，meshVers是原点云，edges为列数为2的边矩阵；

f = fopen( fileName, 'w' );
edgesCount = size(edges, 1);

% 取edges涉及到的顶点；
versIdx = [];
for i = 1: size(edges, 1)
    versIdx = [versIdx, edges(i, :)];
end
chosenVers = meshVers(versIdx, :);

% 生成写入文件中的边数据：[1,2; 3,4; 5,6; ...]
newEdges = [];
for i = 1: edgesCount
    newEdges = [newEdges; [2*i - 1, 2*i]];
end


fprintf( f, 'v %0.17g %0.17g %0.17g\n', chosenVers');
fprintf( f, 'l %d %d \n', newEdges');

end