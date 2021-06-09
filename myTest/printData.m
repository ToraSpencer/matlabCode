clc
clear all
functionname='testforgumline517.m'; functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir 'trianglerayintersection'])
addpath([functiondir 'toothmodel'])
addpath([functiondir 'theratofp'])
addpath([functiondir 'Root deformation'])
addpath([functiondir 'modelread'])
addpath([functiondir 'mixFE'])
addpath([functiondir 'margin line'])
addpath([functiondir 'Arch curve'])
addpath([functiondir 'Dentalmodelwithaxis'])
addpath([functiondir '������루��׼���Ͳ������ڶ��룩'])
addpath([functiondir 'siofmodeltooth'])
addpath([functiondir 'attachments'])
addpath([functiondir 'modelStadard'])
addpath([functiondir 'testdata'])
addpath([functiondir 'mesh process'])
addpath([functiondir 'testdata(4)'])
% % ���ر�׼��
load('dental_crown.mat');
load('dentalmodelwithroot0.1forlow.mat');
load('axisofdental_crowm');
load('trfromtoothtocrown.mat');
 
 
for x = 11:49
 
        
s = fix(x/10);  %ȡ��
g = mod(x,10);%ȡ��


if s == 1 || s == 2
    fdi = textread('FDIUpper__.dxt');
    toothIdx = find (fdi == x);

    if isempty(toothIdx)
        disp( '�ò���û�д����� ');
 
    else
        fdi = textread('FDIUpper__.dxt');
        toothIdx = find (fdi == x);
        namestr2 = ['gumlineUpper_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);
        str = ['gumline', num2str(x), '.obj'];
        OBJwriteVertices(str, gumline);
    end
   
else
    fdi = textread('FDILower__.dxt');
    toothIdx = find (fdi == x);
    if isempty(toothIdx)
        disp( '�ò���û�д����� ');
    else
        fdi = textread('FDILower__.dxt');
        toothIdx = find (fdi == x);
        namestr2 = ['gumlineLower_',num2str(toothIdx-1),'.','obj'];
        gumline =  ReadObj(namestr2);
        str = ['gumline', num2str(x), '.obj'];
        OBJwriteVertices(str, gumline);
    end
    

        
     
end

end
disp('finished');
  