classdef TeachingMachine

    properties
    end
    
    methods(Static)
        function compare(realNumber1,realNumber2)
            if(realNumber1>realNumber2)
                fprintf('�������нϴ����Ϊ��%f\n\n',realNumber1);
            else
                fprintf('�������нϴ����Ϊ��%f\n\n',realNumber2);
            end
        end
        
        %%��ͨ�����ĺ������
        function functionHandle1(realNumber1,realNumber2)           
            f = @findMax;
            max = f(realNumber1,realNumber2);
            disp('���������������һ��������');
            disp('��ͨ�����ĺ��������findMax()����');
            fprintf('�������нϴ����Ϊ��%f\n\n',max)
        end
        
        %%�෽���ĺ������
        function functionHandle2(realNumber1,realNumber2)           
            f = @TeachingMachine.compare;
            disp('�෽���ĺ��������Static�෽����TeachingMachine.compare()');
            f(realNumber1,realNumber2);

        end
        
        %%�����������
        function functionHandle3()
            g = @(x) x.^2+3;
            disp('����������������ټ��ض��塢ʹ�ú�����');
            fprintf('g(3) = %d\n',g(3));
            fprintf('g(18) = %d\n\n',g(18));

        end
        
    
       
    end
    

end

