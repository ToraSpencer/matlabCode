classdef ATM <handle
    %ATM��
    %���ԣ��������е����ƣ��ֽ����(�����³��ֽ���ͻ�)���ж��ڲ��Ƿ������п���flag����ȡ�����п��Ŀ��ţ���ȡ�����п��Ĵ�
    %������1. �ж��ڲ��Ƿ������п���flag��
    %     2. �������п�
    %     3. �³����п�
    %     4. ��ȡ���п�������
    
 
    properties
        bankName
        ATMcash
        cardFlag
        
    end
    
    methods
        function obj = ATM(bankName)
            obj.bankName = bankName;
            obj.ATMcash = cash(100000);             %��ʼʱATM�����д������ֽ�
            disp(strcat('һ��',obj.bankName,'ATM������ʵ������'));
        end
        
        function eatCard(obj,Card)
            if obj.cardFlag == 0
                type = superclasses(Card);
                if type{1} == 'bankCard'
                    obj.cardFlag = 1;
                    disp(strcat(obj.bankName,'��ӭ����'));
                else
                    disp('�����������п���');
                    obj.spit(Card);
                end  
            obj.loginInterface(Card);                %�����¼���档 
            else
                error('�������п��������޷��忨��');
            end
        end
        
        function loginInterface(obj,Card)
            x = input('��ѡ�������');
            switch x
                case 1
                    password = input('���������룺');
                    if password == Card.password
                        obj.mainInterface(Card);            %������������档
                    else
                        error('�������');
                    end
                        
                case 2
                    obj.spit(Card);
                otherwise
                    error('������������û�ж�Ӧѡ�');
            end
            
        end
        
        function mainInterface(obj,Card)
            if class(Card) == 'DepositCard'
                x = input('��ѡ������Ҫ�ķ���\n1.���\n2.ȡ��\n3.��ѯ\n4.�˿�\n');
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
                        error('����û�ж�Ӧ�ķ���');
                end
                        
            else  
                x = input('��ѡ������Ҫ�ķ���\n1.���\n2.ȡ��\n3.��ѯ\n4.���ÿ�����\n5.�˿�\n');
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
                        error('����û�ж�Ӧ�ķ���');
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
            disp('ATM���³���һ�ſ�');
            obj.cardFlag = 0;
        end
        
        
        
        
        
        
        
    end
    
end

