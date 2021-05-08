classdef vector3D < vector2D       %��ά������̳��ڶ�ά������     
    %�ѿ�������ϵ�¡���ά�ռ��е�������
    %
    
    properties
        z
    end
    
    properties(Dependent)
        r              %ò�����Բ�����д��ֻ�ܶ����µ�����
    end
   
    
    methods
        function obj = vector3D(x,y,z)
            %��д���췽�������ﲻ��ͨ���ȼ̳и��෽�������������ʵ����д������ֱ����д��
            obj = obj@vector2D(x,y);
            if nargin==0
               obj.x = 1;
               obj.y = 0;
               obj.z = 0;
            elseif nargin==3
                obj.x = x;
                obj.y = y;
                obj.z = z;
            end
        end
              
       function set.z(obj,val)
            if isnumeric(val)==1
                obj.z = val;
            else
                error('z��������ֵ��');
            end
        end
        
        
        function r = get.r(obj)
             r = sqrt(obj.x^2+obj.y^2+obj.z^2); 
        end
        
        
        function normalize(obj)
            %��д��һ���������ȼ̳и����еĹ�һ����������������ݡ�
            %���룺һ������
            %������÷���ֻ�Ǹı����ж�������ԣ�������û�з���ֵ�ġ�
            norm = obj.r;
            obj.x = obj.x/norm;
            obj.y = obj.y/norm;
            obj.z = obj.z/norm;
        end
        
        function disp(obj)
            %��дdisp���������ﲻ��ͨ���ȼ̳и��෽�������������ʵ����д������ֱ����д��
            str = sprintf('����������Ϊ��\n');                   %ʹ���ַ��������������ϵ������ַ��������ݡ�
            str = [str,sprintf('(%.4f,%.4f,%.4f)\n',obj.x,obj.y,obj.z)];
            str = [str,sprintf('������ģΪ��\n%.4f\n',obj.r)];
            disp(str);
        end
    end
end

