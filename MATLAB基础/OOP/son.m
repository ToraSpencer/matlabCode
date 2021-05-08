classdef son <father

    properties
    end
    
    methods
        function obj = son(name,age)
            obj = obj@father(name,age);            %子类必须继承父类的构造方法。
            disp('一个儿子对象被实例化了');
        end
        
        function smoke(obj)
            smoke@father(obj);
            disp('子类方法：抽烟被调用了。');
        end
        
        function o = father(obj)
            %儿子→爸爸的转换函数
            o = father(obj.name,obj.age);
            disp('我升级成爸爸啦！');
        end
   
    end
end

