%% �������ݽṹ��

clc;
clear all;
close all;

%% �ṹ��struct
tooth.vers = [];                % ֱ������struct���ֶΣ�
tooth.tris = [];
[tooth.vers, tooth.tris] = readOBJ('./data/tooth.obj');