classdef head <handle
    %ͷ�࣬���۱ǿڶ������ϡ�ͷ�����۱ǿڶ�Ҳ�Ͳ��ܹ����ˣ�����ͷ������ĸ��û�еĻ�����Ϊ���ֵ��۱ǿڶ��ǲ��ܵ��������ġ�
    %���ԣ��۱ǿڶ�����
    %���������췽�������ˣ�������������������
    
    properties                         %��Щ���Ե����Ͷ��Ƕ���
        myEyes
        myNose
        myEars
        myMouth
    end
    
    methods
        function obj = head()           %��Ϊ���ĸ������ͷ�������ˣ��䲿������Ҳ�����������
            obj.myEyes = eyes();
            obj.myNose = nose();
            obj.myEars = ears();
            obj.myMouth = mouth();
        end

        
 
        function say(obj)
            disp('i''m a head.');
        end

         function laugh(obj)
            %���룺ͷ����
            %�����
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

