clc;
clear all;
handle1 = vector2D(2,7);
handle2 = handle1;
%%
%�̳�handle����ʵ���������Ķ��������Ǿ�����������ָ����ͬһ��value���󣬵���һ������ķ����ı������value��������ԣ���һ��������������Ҳ��֮�仯��
handle1.x = 5;
disp(handle1.x);
disp(handle2.x);

%%
%handle�����������������
handle1.delete();
fprintf('%d   %d\n',isvalid(handle1),isvalid(handle2));  %handle��ָ�����ݱ��ͷţ�����handle������Ч����

%%
%��Ч�����handle�����������ָ�����ݡ�
handle1 = vector2D(2,7);
handle2 = handle1;
disp(handle1);          
disp(handle2);

%%
%��handle����ʹ��clear������
clear handle1;
fprintf('%d\n',isvalid(handle2)); 
