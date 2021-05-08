classdef weapon <handle
    
    properties
        name
        length
    end
    
    methods
        function obj = weapon(name,length)
            obj.name = name;
            obj.length = length;
            disp(strcat('一个武器对象：',obj.name,'被实例化了。'));
        end
        

    end
end

