function [ versOut, trisOut ] = selectedTris2Mesh( vers, tris, selectedTrisIdx )
% 原网格中有一系列的三角片被选取，索引存在向量selectedTrisIdx中，输出这些三角片确定的新网格；

versCount = size(vers, 1);
trisCount = size(tris, 1);
newTrisCount = numel(selectedTrisIdx);

selectedTris = tris(selectedTrisIdx, :);
versIdx = (1: versCount)';
flag = ismember(versIdx, selectedTris);
newOldIdxInfo = versIdx(flag);
oldNewIdxInfo = -ones(versCount, 1);

newVersCount = numel(newOldIdxInfo);
for i = 1: newVersCount
    index = newOldIdxInfo(i);
    oldNewIdxInfo(index) = i;
end

% 老的顶点索引映射成新的索引；
versOut = vers(newOldIdxInfo, :);
trisVec = reshape(selectedTris', 3*newTrisCount, 1);       % 三角片数据拉成一条列向量；tris = (reshape(trisVec, 3, trisCount))';
newTrisVec = oldNewIdxInfo(trisVec);
trisOut = (reshape(newTrisVec, 3, newTrisCount))';

end

