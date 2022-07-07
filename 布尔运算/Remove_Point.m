function [newVers, newTris, reidx] = Remove_Point(vers, tris)
    versCount = size(vers,1);
    idx = ismember(1:versCount, tris(:)');                % 未被索引点
    tempVers = vers(idx,:);
    [newVers,~,ic] = unique(tempVers, 'rows');         % 重复点
    reidx = zeros(1,versCount);
    reidx(idx) = ic;
    newTris = reidx(tris);
end
 