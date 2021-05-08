classdef head <handle
    %头类，是眼鼻口耳类的组合。头坏了眼鼻口耳也就不能工作了，所以头这个组合母体没有的话，作为部分的眼鼻口耳是不能单独工作的。
    %属性：眼鼻口耳对象。
    %方法：构造方法，高兴，生气，悲哀，哭泣。
    
    properties                         %这些属性的类型都是对象。
        myEyes
        myNose
        myEars
        myMouth
    end
    
    methods
        function obj = head()           %作为组合母体对象的头被构造了，其部件对象也被构造出来。
            obj.myEyes = eyes();
            obj.myNose = nose();
            obj.myEars = ears();
            obj.myMouth = mouth();
        end

        
 
        function say(obj)
            disp('i''m a head.');
        end

         function laugh(obj)
            %输入：头对象。
            %输出：
            obj.myEyes.laugh();
            obj.myNose.laugh();
            obj.myEars.laugh();
            obj.myMouth.laugh();
        end
         
        function getAngry(obj)
            obj.myEyes.getAngry()
            obj.myNose.getAngry()
            obj.myEars.getAngry()
            obj.myMouth.getAngry()
        end
        
        function getSad(obj)
            obj.myEyes.getSad()
            obj.myNose.getSad()
            obj.myEars.getSad()
            obj.myMouth.getSad()
        end
        
        function cry(obj)
            obj.myEyes.cry()
            obj.myNose.cry()
            obj.myEars.cry()
            obj.myMouth.cry()
        end
        
end    
end

