classdef ATM <handle
    %ATM机
    %属性：所属银行的名称，现金对象(可以吐出现金给客户)，判断内部是否有银行卡的flag，读取的银行卡的卡号，读取的银行卡的存款。
    %方法：1. 判断内部是否有银行卡的flag。
    %     2. 吸入银行卡
    %     3. 吐出银行卡
    %     4. 读取银行卡的数据
    
 
    properties
        bankName
        ATMcash
        cardFlag
        
    end
    
    methods
        function obj = ATM(bankName)
            obj.bankName = bankName;
            obj.ATMcash = cash(100000);             %初始时ATM机内有大量的现金。
            disp(strcat('一个',obj.bankName,'ATM机对象被实例化了'));
        end
        
        function eatCard(obj,Card)
            if obj.cardFlag == 0
                type = superclasses(Card);
                if type{1} == 'bankCard'
                    obj.cardFlag = 1;
                    disp(strcat(obj.bankName,'欢迎您！'));
                else
                    disp('请勿插入非银行卡。');
                    obj.spit(Card);
                end  
            obj.loginInterface(Card);                %进入登录界面。 
            else
                error('机器内有卡，现在无法插卡。');
            end
        end
        
        function loginInterface(obj,Card)
            x = input('请选择操作：');
            switch x
                case 1
                    password = input('请输入密码：');
                    if password == Card.password
                        obj.mainInterface(Card);            %进入服务主界面。
                    else
                        error('密码错误');
                    end
                        
                case 2
                    obj.spit(Card);
                otherwise
                    error('错误：输入数字没有对应选项。');
            end
            
        end
        
        function mainInterface(obj,Card)
            if class(Card) == 'DepositCard'
                x = input('请选择您需要的服务：\n1.存款\n2.取款\n3.查询\n4.退卡\n');
                switch x
                    case 1
                        obj.deposit(Card);
                    case 2
                        obj.debit(Card);
                    case 3
                        obj.check(Card);
                    case 4
                        obj.spit(Card);
                    otherwise
                        error('输入没有对应的服务。');
                end
                        
            else  
                x = input('请选择您需要的服务：\n1.存款\n2.取款\n3.查询\n4.信用卡服务\n5.退卡\n');
                switch x
                    case 1
                        obj.deposit(Card);
                    case 2
                        obj.debit(Card);
                    case 3
                        obj.check(Card);
                    case 4
                        obj.upgrade(Card);
                    case 5
                        obj.spit(Card)
                    otherwise
                        error('输入没有对应的服务。');
                end
            end
            
        end
        
        function deposit(obj,Card)
        end
        
        function debit(obj,Card)
        end
        
        function check(obj,Card)
        end
        
        function upgrade(obj,Card)
        end
        
        function spit(obj)
            disp('ATM机吐出了一张卡');
            obj.cardFlag = 0;
        end
        
        
        
        
        
        
        
    end
    
end

