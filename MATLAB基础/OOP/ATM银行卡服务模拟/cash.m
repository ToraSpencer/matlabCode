classdef cash <handle
    %�ֽ��ࡣ
    %���û����ATM������Ϲ�ϵ�������ǵĲ����ࡣ�ֽ����ܽ����κλ����������û�����塣
    %���ԣ����
    %�����������������ֽ�������Ϊ�����Ի٣�
    
    properties
        Amount
    end
    
    methods
        function obj = cash(Amount)
            obj.Amount = Amount;
        end
        
        function AmountUp(obj,value)
        end
        
        function AmountDown(obj,value)
        end
        
    end
end

