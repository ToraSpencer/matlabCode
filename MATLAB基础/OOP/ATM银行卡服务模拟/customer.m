classdef customer <handle
    %ʹ�����п�����
    %�����п��ࡢ�ֽ�������Ϲ�ϵ���û���ĸ���࣬���Բ������п����ֽ𣬼����Ե������ǵķ�����
    %����:���п������ֽ����
    %������
    %         1.�������п��ķ�����
    %                        1.1.�������п�
    %                        1.2.ȡ���³������п�
    %         2.����ATM�ķ�����
    %                        2.1.�������п�����
    %                        2.2.���ATM�����ϵ�ѡ����1,2,3,..8
    %         3.��ATM��������ֽ�ķ�����
    %         4.ȡ��ATM���³��ֽ�ķ�����
                            
  
    
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
            disp(strcat('ʵ������һ���ͻ�����',obj.name));
        end
        
        function chooseAcard(obj,ATMmachine)       %ѡ������ſ���
            x = input('��ѡ��Ҫ�������ſ���\n1.���\n2.���ÿ�\n');
            if x == 1
                ATMmachine.eatCard(obj.myDepositCard);
            elseif x==2
                ATMmachine.eatCard(obj.myCreditCard);
            else
                error('������ѡ��1��2��');
            end
        end
        
        function putInCash(obj,ATMmachine,value)         %��ATM��������ֽ�
            obj.myCash.AmountDown(value);
            disp(strcat(obj.name,'�������ֽ𣬵����ȷ�ϡ�'));
            
            
            
        end
        
        
            
            
            
            
   
        

        
        
    end
    
end

