classdef bankCard <handle
    %银行卡：分为信用卡和储蓄卡，信用卡可以透支，储蓄卡不可以透支。
    %是一个抽象类。是信用卡类和储蓄卡类的父类。
    %和用户类是组合关系，是用户类的部件，没有用户银行卡本身做不了任何事情，没有存在的意义。
    %属性：表明自己是插入ATM状态还是离开ATM状态的FLAG，卡号，密码，存款数额，额度。
    %方法：插入ATM，离开ATM，存款数额的加减。
    
    properties
        flag
        account
        password
        deposit
        limite
    end
    
    methods
        function obj = bankCard(deposit)
            obj.deposit = deposit;
        end
        
        function insert(obj)        %插卡操作。
            obj.flag = 1;
        end

        function remove(obj)        %拔卡操作。
            obj.flag = 0;
        end
        
        function depositMoney(obj,depositAmount)    %存钱操作。
            obj.deposit = obj.deposit+depositAmount;
        end
        
        function success = debitMoney(obj,debitAmount)  %取钱操作，操作成功返回1
            if obj.deposit-debitAmount>=obj.limit
                obj.deposit = obj.deposit-debitAmount;
                success = 1;
            else
                success = 0;
            end
            
        end
        
        

        
        

        
    end
    
end

