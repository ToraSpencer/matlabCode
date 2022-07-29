%% ��ݼ�
% �۵�����    CTRL + =   
% ȡ���۵�    CTRL + SHIFT + =
% ����ע��    CTRL + R
% ���з�ע��  CTRL + T
% �򿪱��������� CTRL + D
%           �򿪺���ָ���Ǵ򿪺����ļ���
%           ����ֻ����debug��ʱ��򿪣�������ʾ���ݵ����顣


%% �����㷨����sort()
clc;
clear all;
x = round(100 * rand(1, 10));
disp(x);
[u, v] = sort(x);
disp(u);
disp(v);


%% ��ͼ
clc;
clear all;
close all;
x = 1:100;
y = 0.1*(x-35).^2 + 10;

figure('Name', 'xy graph');     

% plot()����������ͼ
plot(x, y, 'b-*');
xlabel('axis X');
ylabel('axis Y');
hold on                  % ������ǰ�������еĻ�ͼ���Ӷ�ʹ����ӵ��������еĻ�ͼ����ɾ�����л�ͼ;

% scatter() ������ɢ��ͼ
x1 = 10*[1: 10];
y1 = -0.1*(x1-40).^2 + 30;
y2 = y1-30;
scatter(x1, y1, 'r');                       % ���ĵ�Ĭ���ǿ���Բ��
scatter(x1, y2, 'g', 'filled');             % filledָ�����ĵ�Ϊʵ��Բ��

% surf()���������棺
% meshgrid()��������դ��
% xlabel(), ylabel�����趨����Ϣ
[X,Y]=meshgrid(-2:0.1:2, -3: 0.1 :3);
Z=exp(-X.^2-Y.^2);
figure(2);
xlabel('x');
ylabel('y');
%           �޸�ͼ���������ԣ�
surfHandle = surf(X,Y,Z);          
surfHandle.EdgeColor = 'none';              % ���رߣ�
        


%% ��д.dat�ļ�
clc;
clear all;
load('./data/elliSampleVers.mat');
save ./data/elliSampleVers2.mat sampleVers;
x = sampleVers(:, 1);
y = sampleVers(:, 2);
figure('Name', 'xy graph');
plot(x, y, 'b-*');


%% ���ÿ���̨����
clc;
clear all;
fileName = 'G:/gitRepositories/matlabCode/MATLAB����/data/tooth.obj';
tic
command = ['E:/workstation/smarteeproj/Exe/Win32_Release/planeCutMesh.exe', ...
   ' planeCutMesh ', fileName, ' [0, 0, 0] [0, 0, 1]'];
[status, result] = system( command );
fprintf('planeCutMesh calculation takes %f s time.\n', toc);


 