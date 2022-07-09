function objWriteEdges( fileName, edges, meshVers)
% �����ʾ�ߵ�OBJ�ļ���meshVers��ԭ���ƣ�edgesΪ����Ϊ2�ı߾���

f = fopen( fileName, 'w' );
edgesCount = size(edges, 1);

% ȡedges�漰���Ķ��㣻
versIdx = [];
for i = 1: size(edges, 1)
    versIdx = [versIdx, edges(i, :)];
end
chosenVers = meshVers(versIdx, :);

% ����д���ļ��еı����ݣ�[1,2; 3,4; 5,6; ...]
newEdges = [];
for i = 1: edgesCount
    newEdges = [newEdges; [2*i - 1, 2*i]];
end


fprintf( f, 'v %0.17g %0.17g %0.17g\n', chosenVers');
fprintf( f, 'l %d %d \n', newEdges');

end