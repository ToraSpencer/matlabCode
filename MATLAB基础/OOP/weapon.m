classdef weapon <handle
    
    properties
        name
        length
    end
    
    methods
        function obj = weapon(name,length)
            obj.name = name;
            obj.length = length;
            disp(strcat('һ����������',obj.name,'��ʵ�����ˡ�'));
        end
        

    end
end

