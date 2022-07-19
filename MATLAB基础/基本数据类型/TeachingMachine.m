classdef TeachingMachine

    properties
    end
    
    methods(Static)
        function compare(realNumber1,realNumber2)
            if(realNumber1>realNumber2)
                fprintf('两个数中较大的数为：%f\n\n',realNumber1);
            else
                fprintf('两个数中较大的数为：%f\n\n',realNumber2);
            end
        end
        
        %%普通函数的函数句柄
        function functionHandle1(realNumber1,realNumber2)           
            f = @findMax;
            max = f(realNumber1,realNumber2);
            disp('函数句柄：给函数一个别名。');
            disp('普通函数的函数句柄：findMax()函数');
            fprintf('两个数中较大的数为：%f\n\n',max)
        end
        
        %%类方法的函数句柄
        function functionHandle2(realNumber1,realNumber2)           
            f = @TeachingMachine.compare;
            disp('类方法的函数句柄：Static类方法：TeachingMachine.compare()');
            f(realNumber1,realNumber2);

        end
        
        %%匿名函数句柄
        function functionHandle3()
            g = @(x) x.^2+3;
            disp('匿名函数句柄：快速简介地定义、使用函数。');
            fprintf('g(3) = %d\n',g(3));
            fprintf('g(18) = %d\n\n',g(18));

        end
        
    
       
    end
    

end

