classdef customer <handle
    %使用银行卡的人
    %和银行卡类、现金类是组合关系。用户是母体类，可以操作银行卡和现金，即可以调用它们的方法。
    %属性:银行卡对象，现金对象。
    %方法：
    %         1.操作银行卡的方法：
    %                        1.1.插入银行卡
    %                        1.2.取走吐出的银行卡
    %         2.操作ATM的方法：
    %                        2.1.输入银行卡密码
    %                        2.2.点击ATM界面上的选项：编号1,2,3,..8
    %         3.往ATM机里放入现金的方法。
    %         4.取走ATM机吐出现金的方法。
                            
  
    
    properties
        name
        myDepositCard
        myCreditCard
        myCash
    end
    
    methods
        function obj = customer(name,originalDeposit1,originalDeposit2,originalCash)
            obj.name = name;
            obj.myDepositCard = DepositCard(originalDeposit1);
            obj.myCreditCard = CreditCard(originalDeposit2);
            obj.myCash = cash(originalCash);
            disp(strcat('实例化了一个客户对象：',obj.name));
        end
        
        function chooseAcard(obj,ATMmachine)       %选择插哪张卡。
            x = input('请选择要插入哪张卡：\n1.储蓄卡\n2.信用卡\n');
            if x == 1
                ATMmachine.eatCard(obj.myDepositCard);
            elseif x==2
                ATMmachine.eatCard(obj.myCreditCard);
            else
                error('错误：请选择1或2。');
            end
        end
        
        function putInCash(obj,ATMmachine,value)         %往ATM机里放入现金。
            obj.myCash.AmountDown(value);
            disp(strcat(obj.name,'放入了现金，点击了确认。'));
            
            
            
        end
        
        
            
            
            
            
   
        

        
        
    end
    
end

