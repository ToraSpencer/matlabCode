%% ��ȡ����ָʾ�ߡ�������ָʾ����ĵ�λ�������������ָʾ�ߵ㼯

function [ dirLine ] = getDirLine( dir, start )
dir = dir./norm(dir, 2);
SR = 0.5;			% �ռ������SR������������������ľ��루��λmm��
length = 10;
versCount = round(length / SR);
dirLine = ones(versCount, 3);
dirLine(1, :) = start;
for k = 2 : versCount
    dirLine(k, :) = dirLine(1, :) + SR * dir * (k-1);
end



end

