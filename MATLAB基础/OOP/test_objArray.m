clc;
clear all;

%%
%�����������鷽��1��
personArray1(1) = person('a',2);
personArray1(2) = person('b',1);
personArray1(3) = person('c',23);

%%
%�����������鷽��2��
person1 = person('d',12);
person2 = person('e',23);
person3 = person('f',12);
personArray2 = [person1,person2,person3]; 



%%
%��ȡ�������顣
personArray3 = personArray2(1:2);

%%
%�������飺��������û���Լ���delete()������ֻ��������Ԫ�ض������delete���������Լ���
disp(personArray1(1));
for i = 1:length(personArray1)
    personArray1(i).delete();
end
disp(personArray1(1));      %handle�������ٺ���Ȼ���ڣ�ֻ�Ǳ���Ч����
disp(personArray1);         %handle�������������Ԫ�ض��������󣬶���������Ȼ����

%%
%
personArray1 = personArray2;








