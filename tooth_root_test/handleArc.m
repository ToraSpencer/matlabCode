function [ arcOut ] = handleArc( arcIn )
% �����⻡��ת��Ϊ[-pi, pi)�Ļ���
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

