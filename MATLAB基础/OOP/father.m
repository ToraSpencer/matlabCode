classdef father <person
    

    properties
        myfoot
    end
    
    methods
        function obj = father(name,age)
             obj = obj@person(name,age);
            obj.myfoot = foot();
            disp('һ�����׶���ʵ�����ˡ�');
        end
        
        function smoke(obj)
            disp('���෽�������̱������ˡ�');
        end
        
        
        function fatherBeatSon1(obj,skillName)       %ĸ�������ò����������Ժͷ�����
            obj.myfoot.uniqueSkill = skillName;
            disp(strcat('����ѡ���˱�ɱ����',skillName));
            disp(strcat('����������',obj.myfoot.uniqueSkill,'�����'));
        end
               
        function fatherBeatSon2(obj,weapon)          %�������Ϊһ������ķ�����
            disp(strcat('���������ó�',num2str(weapon.length),'�׵�',weapon.name,'����ӡ�'));
        end
        
        function smokeAndBeat(obj)
            obj.smoke();                             %������Ƕ�ס�
            disp('����һ�߳���һ�ߴ���ӡ�');
        end
        
        
        
        
        
    end
end

