classdef vector2D <handle      %<handle一定要写！！！
    %笛卡尔坐标系下，二维平面上的向量。
    % 
    
    properties
        x                %只声明属性，不设默认值。
        y 
        objName
    end
    
    properties(Dependent)    %非独立属性，有该属性时，一定要写get方法。
        rou                  %声明就可以了，不需要初始化。
    end

    methods
        function obj = vector2D(x,y)
            %构造方法：将类实例化输出对象的方法
            %输入：想要实例化的对象的独立属性。也可以没有输入，那就是default constructor。
            %输出：一个对象。
            if nargin==2
                obj.x = x;
                obj.y = y;
            else
                obj.x = 1;
                obj.y = 1;
            end  
        end
        
        function rou = get.rou(obj)
            %get方法：访问对象的属性。
            %输入：一个对象。
            %输出：一个属性的值。
             rou = sqrt(obj.x^2+obj.y^2); 
        end 
        
        function normalize(obj)
            %归一化方法。
            %输入：一个对象。
            %输出：无输出。该方法只是改变已有对象的属性，所以是没有返回值的。
            norm = obj.rou;          %不能把obj.rou直接带入到下面的计算式中，因为obj.x一旦改变，obj.rou也就立刻随之改变
            obj.x = obj.x/norm;
            obj.y = obj.y/norm;
        end
      
        function disp(obj)
            %重写disp方法，使得对象以我们想要的形式输出。
            str = sprintf('向量的坐标为：\n');                   %使用字符向量串联，不断地增加字符串的内容。
            str = [str,sprintf('(%.4f,%.4f)\n',obj.x,obj.y)];
            str = [str,sprintf('向量的模为：\n%.4f\n',obj.rou)];
            disp(str);
       end
    end
    
    
methods(Sealed,Static)              %封口的方法，无法被继承。
     function say()
           disp('i''m a vector!');
       end
end


    
end

