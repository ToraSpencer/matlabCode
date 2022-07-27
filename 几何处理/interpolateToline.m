function [ line ] = interpolateToline( startVer, endVer, SR)
% 两点生成一条直线点集
 dir = endVer - startVer;
 length = norm(dir);
 if(length < SR)
     line = [];
 else
     dir = normalizerow(dir);
     offsetVecs = [0 : SR : length]' * dir;
     segCount = size(offsetVecs, 1);
     starts = repmat(startVer, segCount, 1);
     line = starts + offsetVecs;
     dis = endVer - line(segCount, :);
     
     % 插入的最后一个点和endVers距离不可以太近；
     if(norm(dis) < 0.8 * SR)
         line(segCount, :) = endVer;
     else
         line = [line; endVer];
     end
     

 end
end

