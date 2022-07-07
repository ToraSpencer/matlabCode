%% 获取方向指示线——输入指示方向的单位向量，输出方向指示线点集

function [ dirLine ] = getDirLine( dir, start )
dir = dir./norm(dir, 2);
SR = 0.5;			% 空间采样率SR——相邻两个采样点的距离（单位mm）
length = 10;
versCount = round(length / SR);
dirLine = ones(versCount, 3);
dirLine(1, :) = start;
for k = 2 : versCount
    dirLine(k, :) = dirLine(1, :) + SR * dir * (k-1);
end



end

