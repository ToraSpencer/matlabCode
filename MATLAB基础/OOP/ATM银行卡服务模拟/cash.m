classdef cash <handle
    %现金类。
    %和用户类和ATM类是组合关系，是它们的部件类。现金本身不能进行任何活动，单独存在没有意义。
    %属性：金额
    %方法：析构函数（现金用完视为对象自毁）
    
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

