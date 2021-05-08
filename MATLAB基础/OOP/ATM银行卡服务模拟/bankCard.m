classdef bankCard <handle
    %���п�����Ϊ���ÿ��ʹ�������ÿ�����͸֧�����������͸֧��
    %��һ�������ࡣ�����ÿ���ʹ����ĸ��ࡣ
    %���û�������Ϲ�ϵ�����û���Ĳ�����û���û����п������������κ����飬û�д��ڵ����塣
    %���ԣ������Լ��ǲ���ATM״̬�����뿪ATM״̬��FLAG�����ţ����룬��������ȡ�
    %����������ATM���뿪ATM���������ļӼ���
    
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
        
        function insert(obj)        %�忨������
            obj.flag = 1;
        end

        function remove(obj)        %�ο�������
            obj.flag = 0;
        end
        
        function depositMoney(obj,depositAmount)    %��Ǯ������
            obj.deposit = obj.deposit+depositAmount;
        end
        
        function success = debitMoney(obj,debitAmount)  %ȡǮ�����������ɹ�����1
            if obj.deposit-debitAmount>=obj.limit
                obj.deposit = obj.deposit-debitAmount;
                success = 1;
            else
                success = 0;
            end
            
        end
        
        

        
        

        
    end
    
end

