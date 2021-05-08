classdef father <person
    

    properties
        myfoot
    end
    
    methods
        function obj = father(name,age)
             obj = obj@person(name,age);
            obj.myfoot = foot();
            disp('一个父亲对象被实例化了。');
        end
        
        function smoke(obj)
            disp('父类方法：抽烟被调用了。');
        end
        
        
        function fatherBeatSon1(obj,skillName)       %母类对象调用部件对象属性和方法。
            obj.myfoot.uniqueSkill = skillName;
            disp(strcat('父亲选择了必杀技：',skillName));
            disp(strcat('父亲正在用',obj.myfoot.uniqueSkill,'打儿子'));
        end
               
        function fatherBeatSon2(obj,weapon)          %输入参数为一个对象的方法。
            disp(strcat('父亲正在用长',num2str(weapon.length),'米的',weapon.name,'打儿子。'));
        end
        
        function smokeAndBeat(obj)
            obj.smoke();                             %方法的嵌套。
            disp('父亲一边抽烟一边打儿子。');
        end
        
        
        
        
        
    end
end

