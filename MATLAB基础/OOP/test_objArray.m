clc;
clear all;

%%
%创建对象数组方法1：
personArray1(1) = person('a',2);
personArray1(2) = person('b',1);
personArray1(3) = person('c',23);

%%
%创建对象数组方法2：
person1 = person('d',12);
person2 = person('e',23);
person3 = person('f',12);
personArray2 = [person1,person2,person3]; 



%%
%截取对象数组。
personArray3 = personArray2(1:2);

%%
%销毁数组：对象数组没有自己的delete()方法，只能依次让元素对象调用delete方法析构自己。
disp(personArray1(1));
for i = 1:length(personArray1)
    personArray1(i).delete();
end
disp(personArray1(1));      %handle对象销毁后依然存在，只是被无效化。
disp(personArray1);         %handle对象数组的所有元素对象析构后，对象数组依然存在

%%
%
personArray1 = personArray2;








