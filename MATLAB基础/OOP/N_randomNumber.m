classdef N_randomNumber < handle
    %正态分布随机数对象：normal distributed random number
    properties
        miu = 0
        cumulativeProbability
        value1       
    end
    
    properties(Dependent)
        value2
    end
    
    properties(Constant)
        sigma = 1
        variance = N_randomNumber.sigma^2
    end
    
    methods
        function obj = N_randomNumber()
            obj.value1 = randn();        
        end
        
        function value2 = get.value2(obj)
            value2 = normrnd(obj.miu,obj.sigma);
        end
    end
end

