%���������

%%
%�Զ��庯���ĺ��������   
%               �������f = @��������
%               �������f = str2func('������');
clc;
close all;
clear all;

f = @findMax;
max = f(3,4);
disp('���������������һ��������');
fprintf('�������нϴ����Ϊ��%f\n\n',max)

g = str2func('findMax');
max2 = g(5,6);
fprintf('�������нϴ����Ϊ��%f\n\n',max2)



%%
%������������� �������f = @(�����б�)�������ʽ��
clc;
close all;
clear all;

g = @(x) x.^2+3;
disp('����������������ټ��ض��塢ʹ�ú�����');
fprintf('g(3) = %d\n',g(3));
fprintf('g(18) = %d\n\n',g(18));

%%
%���������Ϊ��һ�������Ĳ�����
clc;
clear all;

ezplot(@sin, [0, 2 * pi]);      %�������@sin��Ϊ���������˺���ezplot;




%%
%�ٶȲ��ԣ�
clc;
close all;
clear all;

M = 1E6;

tic;
y = sinfun1(M);
t1 = toc;
disp('sinfun1��û��Ԥ�ȿ����ڴ��forѭ��');
fprintf('sinfun1()ֱ������ʱ�䣺t1 = %.3g\n',t1);

f1 = @() sinfun1(M);            %M�ǳ��������������������б�Ϊ�ա�
t11 = timeit(f1);               %timeit(�������)���غ���ִ��ʹ��ʱ�䡣
fprintf('sinfun1()ʹ�ú����������ʱ�䣺t11 = %.3g\n',t11);
                                %���ú�������ȵ��ú����ٶȿ졣

tic;
y = sinfun2(M);
t2 = toc;
disp('sinfun2��Ԥ�ȿ����ڴ��˵�forѭ��');              %Ԥ�ȿ����ڴ�ռ����ʱ�����ٶȿ졣
fprintf('sinfun2()ֱ������ʱ�䣺t2 = %.3g\n',t2);

f2 = @() sinfun2(M);            
t21 = timeit(f2);
disp('sinfun3�ǽ�forѭ�������������˵ĺ���')             %����������forѭ������������
fprintf('sinfun2()ʹ�ú����������ʱ�䣺t21 = %.3g\n',t21);

tic;
y = sinfun3(M);                             
t3 = toc;
fprintf('sinfun3()ֱ������ʱ�䣺t3 = %.3g\n',t3);

f3 = @() sinfun3(M);            
t31 = timeit(f3);
fprintf('sinfun3()ʹ�ú����������ʱ�䣺t31 = %.3g\n',t31);