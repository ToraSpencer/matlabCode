classdef vector3D < vector2D       %三维向量类继承于二维向量类     
    %笛卡尔坐标系下、三维空间中的向量。
    %
    
    properties
        z
    end
    
    properties(Dependent)
        r              %貌似属性不能重写，只能定义新的属性
    end
   
    
    methods
        function obj = vector3D(x,y,z)
            %重写构造方法，这里不能通过先继承父类方法再添加内容来实现重写，所以直接重写。
            obj = obj@vector2D(x,y);
            if nargin==0
               obj.x = 1;
               obj.y = 0;
               obj.z = 0;
            elseif nargin==3
                obj.x = x;
                obj.y = y;
                obj.z = z;
            end
        end
              
       function set.z(obj,val)
            if isnumeric(val)==1
                obj.z = val;
            else
                error('z必须是数值。');
            end
        end
        
        
        function r = get.r(obj)
             r = sqrt(obj.x^2+obj.y^2+obj.z^2); 
        end
        
        
        function normalize(obj)
            %重写归一化方法：先继承父类中的归一化方法，再添加内容。
            %输入：一个对象。
            %输出：该方法只是改变已有对象的属性，所以是没有返回值的。
            norm = obj.r;
            obj.x = obj.x/norm;
            obj.y = obj.y/norm;
            obj.z = obj.z/norm;
        end
        
        function disp(obj)
            %重写disp方法，这里不能通过先继承父类方法再添加内容来实现重写，所以直接重写。
            str = sprintf('向量的坐标为：\n');                   %使用字符向量串联，不断地增加字符串的内容。
            str = [str,sprintf('(%.4f,%.4f,%.4f)\n',obj.x,obj.y,obj.z)];
            str = [str,sprintf('向量的模为：\n%.4f\n',obj.r)];
            disp(str);
        end
    end
end

