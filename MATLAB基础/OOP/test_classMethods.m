clc;
clear all;
%%
father = father('�ְ�',40);
son = son('����',12);
weapon1 = weapon('���컭�',30);
weapon2 = weapon('̫�����ڽ�',50);

%%
father.fatherBeatSon1('��������jio'); %ĸ�������ò����������Ժͷ�����

%%
father.fatherBeatSon2(weapon1);     %�������Ϊһ������ķ�����

%%
father.smokeAndBeat();              %�����п���Ƕ�׷�����

%%
son.smoke();                        %����̳и���ķ�����

%%
son.isAlive();                      %�������Static������
person.isAlive();                   %��ֱ�ӵ���Static������

%%
son = son.father();                 %����ת��������ת����������͡�
disp(class(son));


