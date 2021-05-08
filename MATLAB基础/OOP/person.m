classdef person <handle
    properties
        age
        name
    end
    
    properties(Constant)
         livingStatus = 1
    end
    
    
    
    methods
        function obj = person(name,age)
            obj.name = name;
            obj.age = age;
        end
        
        function introduce(obj)
            fprintf(strcat('����',obj.name));
        end
        
        function foo(obj)
            obj.age = 1;
        end
    end
    
    
    
    methods(Static)
       function isAlive()
            if person.livingStatus == 1
                disp('�һ���');
            end
       end
       
       
       
    end
    
    
end

