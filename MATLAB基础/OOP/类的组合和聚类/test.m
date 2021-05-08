%？？？需要学习句柄相关的知识，了解母体对象如何调用部件对象的方法。
%*.尝试销毁母体头对象，看看部件对象还是否存在。
clc;
close all;
clear all;

head = head();
disp(head);

%%
%给我笑！
head.laugh();

%%
%给我生气
head.getAngry();

%%
%给我难受
head.getSad();
%%
%给我哭
head.cry();
