clc;
clear all;
%%
father = father('爸爸',40);
son = son('儿子',12);
weapon1 = weapon('方天画戟',30);
weapon2 = weapon('太乙破阙剑',50);

%%
father.fatherBeatSon1('夺命剪刀jio'); %母类对象调用部件对象属性和方法。

%%
father.fatherBeatSon2(weapon1);     %输入参数为一个对象的方法。

%%
father.smokeAndBeat();              %方法中可以嵌套方法。

%%
son.smoke();                        %子类继承父类的方法。

%%
son.isAlive();                      %对象调用Static方法。
person.isAlive();                   %类直接调用Static方法。

%%
son = son.father();                 %调用转换方法来转换对象的类型。
disp(class(son));


