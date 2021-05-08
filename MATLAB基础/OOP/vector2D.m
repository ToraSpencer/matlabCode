classdef vector2D <handle      %<handleһ��Ҫд������
    %�ѿ�������ϵ�£���άƽ���ϵ�������
    % 
    
    properties
        x                %ֻ�������ԣ�����Ĭ��ֵ��
        y 
        objName
    end
    
    properties(Dependent)    %�Ƕ������ԣ��и�����ʱ��һ��Ҫдget������
        rou                  %�����Ϳ����ˣ�����Ҫ��ʼ����
    end

    methods
        function obj = vector2D(x,y)
            %���췽��������ʵ�����������ķ���
            %���룺��Ҫʵ�����Ķ���Ķ������ԡ�Ҳ����û�����룬�Ǿ���default constructor��
            %�����һ������
            if nargin==2
                obj.x = x;
                obj.y = y;
            else
                obj.x = 1;
                obj.y = 1;
            end  
        end
        
        function rou = get.rou(obj)
            %get���������ʶ�������ԡ�
            %���룺һ������
            %�����һ�����Ե�ֵ��
             rou = sqrt(obj.x^2+obj.y^2); 
        end 
        
        function normalize(obj)
            %��һ��������
            %���룺һ������
            %�������������÷���ֻ�Ǹı����ж�������ԣ�������û�з���ֵ�ġ�
            norm = obj.rou;          %���ܰ�obj.rouֱ�Ӵ��뵽����ļ���ʽ�У���Ϊobj.xһ���ı䣬obj.rouҲ��������֮�ı�
            obj.x = obj.x/norm;
            obj.y = obj.y/norm;
        end
      
        function disp(obj)
            %��дdisp������ʹ�ö�����������Ҫ����ʽ�����
            str = sprintf('����������Ϊ��\n');                   %ʹ���ַ��������������ϵ������ַ��������ݡ�
            str = [str,sprintf('(%.4f,%.4f)\n',obj.x,obj.y)];
            str = [str,sprintf('������ģΪ��\n%.4f\n',obj.rou)];
            disp(str);
       end
    end
    
    
methods(Sealed,Static)              %��ڵķ������޷����̳С�
     function say()
           disp('i''m a vector!');
       end
end


    
end

