%% 基本数据结构：

clc;
clear all;
close all;

%% 结构体struct
tooth.vers = [];                % 直接声明struct的字段；
tooth.tris = [];
[tooth.vers, tooth.tris] = readOBJ('./data/tooth.obj');