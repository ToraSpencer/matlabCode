function [ versOut, trisOut ] = selectedTris2Mesh( vers, tris, selectedTrisIdx )
% ԭ��������һϵ�е�����Ƭ��ѡȡ��������������selectedTrisIdx�У������Щ����Ƭȷ����������

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

% �ϵĶ�������ӳ����µ�������
versOut = vers(newOldIdxInfo, :);
trisVec = reshape(selectedTris', 3*newTrisCount, 1);       % ����Ƭ��������һ����������tris = (reshape(trisVec, 3, trisCount))';
newTrisVec = oldNewIdxInfo(trisVec);
trisOut = (reshape(newTrisVec, 3, newTrisCount))';

end

