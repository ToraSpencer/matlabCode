function [ line ] = interpolateToline( startVer, endVer, SR)
% ��������һ��ֱ�ߵ㼯
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
     
     % ��������һ�����endVers���벻����̫����
     if(norm(dis) < 0.8 * SR)
         line(segCount, :) = endVer;
     else
         line = [line; endVer];
     end
     

 end
end

