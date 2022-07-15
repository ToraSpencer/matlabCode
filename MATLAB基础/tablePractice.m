%用构造函数table()来创建表对象。
clc;
clear all;
name = {'alice';'bob';'cindy'};
scores = [88;89;100];
colName = {'name','scores'};

reportCard = table(name,scores,'VariableNames',colName);
disp(reportCard);

%%
%用转换函数array2table来生成表对象。
clc;
clear all;
scores = 100*rand(5,3);
colName = {'Chinese','Math','English'};

reportCard = array2table(scores,'VariableNames',colName);
disp(reportCard);

%%
%使用"点+变量名"的语法来访问表对象的属性（比如取某一列的数据），返回矩阵或者元胞
temp = reportCard.Chinese;
temp1 = reportCard.Properties.VariableNames;
disp('语文成绩表为：');
disp(temp);
disp(class(temp));
disp('表头为：');
disp(temp1);
disp(class(temp1));


%%
%使用圆括号来访问表对象的数据，返回table
temp = reportCard(1,:);
disp('alice的成绩为：');
disp(temp);
disp(class(temp));

%%
%使用花括号来访问表对象的数据，返回矩阵或元胞。
temp = reportCard{1,:};
disp('alice的成绩为');
disp(temp);
disp(class(temp));


%%
%使用sortrows函数，参考表中的某一列数据，来对表进行排序
clc;
clear all;
name = {'Alice';'Bob';'Cindy';'David';'Eve';'Frank'};
Chinese = [14;14;99;99;99;67];
Math = [14;55;99;12;35;67];
English = [14;55;99;12;35;67];
colName = {'name','Chinese','Math','English'};

reportCard = table(name,Chinese,Math,English,'VariableNames',colName);

%最看重语文成绩，按照语文成绩降序排列。其次看中数学成绩，在语文成绩相同的情况下，按数学成绩降序排列。
reportCard = sortrows(reportCard,{'Chinese','Math'},{'descend','descend'});
disp(reportCard);






