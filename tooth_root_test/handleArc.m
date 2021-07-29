function [ arcOut ] = handleArc( arcIn )
% 将任意弧度转换为[-pi, pi)的弧度
arc = arcIn;

while(arc < -pi || arc >= pi)
    
    if arc >= pi
        arc = arc - 2*pi;
    end
    
    if arc < -pi
        arc = arc + 2*pi;
    end
    
end

arcOut = arc;

end

