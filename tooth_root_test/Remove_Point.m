function [Vnew,Fnew,reidx] = Remove_Point(V,F)
    nV = size(V,1);
    idx = ismember(1:nV, F(:)'); % 未被索引点
    Vtmp = V(idx,:);
    [Vnew,~,ic] = unique(Vtmp, 'rows'); % 重复点
    reidx = zeros(1,nV);
    reidx(idx) = ic;
    Fnew = reidx(F);
end

% function [NewPoint,NewFace] = Remove_Point(Point,Face)
%     nP = size(Point,1);
%     RemainInd = ismember((1:nP)', Face(:));
%     NewPoint = Point(RemainInd,:);
%     
%     NewInd = zeros(nP,1);
%     NewInd(RemainInd) = (1:size(NewPoint,1))';
%     NewFace = NewInd(Face);
% end
