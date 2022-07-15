%�ù��캯��table()�����������
clc;
clear all;
name = {'alice';'bob';'cindy'};
scores = [88;89;100];
colName = {'name','scores'};

reportCard = table(name,scores,'VariableNames',colName);
disp(reportCard);

%%
%��ת������array2table�����ɱ����
clc;
clear all;
scores = 100*rand(5,3);
colName = {'Chinese','Math','English'};

reportCard = array2table(scores,'VariableNames',colName);
disp(reportCard);

%%
%ʹ��"��+������"���﷨�����ʱ��������ԣ�����ȡĳһ�е����ݣ������ؾ������Ԫ��
temp = reportCard.Chinese;
temp1 = reportCard.Properties.VariableNames;
disp('���ĳɼ���Ϊ��');
disp(temp);
disp(class(temp));
disp('��ͷΪ��');
disp(temp1);
disp(class(temp1));


%%
%ʹ��Բ���������ʱ��������ݣ�����table
temp = reportCard(1,:);
disp('alice�ĳɼ�Ϊ��');
disp(temp);
disp(class(temp));

%%
%ʹ�û����������ʱ��������ݣ����ؾ����Ԫ����
temp = reportCard{1,:};
disp('alice�ĳɼ�Ϊ');
disp(temp);
disp(class(temp));


%%
%ʹ��sortrows�������ο����е�ĳһ�����ݣ����Ա��������
clc;
clear all;
name = {'Alice';'Bob';'Cindy';'David';'Eve';'Frank'};
Chinese = [14;14;99;99;99;67];
Math = [14;55;99;12;35;67];
English = [14;55;99;12;35;67];
colName = {'name','Chinese','Math','English'};

reportCard = table(name,Chinese,Math,English,'VariableNames',colName);

%������ĳɼ����������ĳɼ��������С���ο�����ѧ�ɼ��������ĳɼ���ͬ������£�����ѧ�ɼ��������С�
reportCard = sortrows(reportCard,{'Chinese','Math'},{'descend','descend'});
disp(reportCard);






