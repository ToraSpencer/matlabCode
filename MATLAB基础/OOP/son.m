classdef son <father

    properties
    end
    
    methods
        function obj = son(name,age)
            obj = obj@father(name,age);            %�������̳и���Ĺ��췽����
            disp('һ�����Ӷ���ʵ������');
        end
        
        function smoke(obj)
            smoke@father(obj);
            disp('���෽�������̱������ˡ�');
        end
        
        function o = father(obj)
            %���ӡ��ְֵ�ת������
            o = father(obj.name,obj.age);
            disp('�������ɰְ�����');
        end
   
    end
end

