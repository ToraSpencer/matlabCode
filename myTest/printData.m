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
addpath([functiondir '牙轴对齐（标准牙和病人牙冠对齐）'])
addpath([functiondir 'siofmodeltooth'])
addpath([functiondir 'attachments'])
addpath([functiondir 'modelStadard'])
addpath([functiondir 'testdata'])
addpath([functiondir 'mesh process'])
addpath([functiondir 'testdata(4)'])
% % 加载标准牙
load('dental_crown.mat');
load('dentalmodelwithroot0.1forlow.mat');
load('axisofdental_crowm');
load('trfromtoothtocrown.mat');
 
 
for x = 11:49
 
        
s = fix(x/10);  %取整
g = mod(x,10);%取余


if s == 1 || s == 2
    fdi = textread('FDIUpper__.dxt');
    toothIdx = find (fdi == x);

    if isempty(toothIdx)
        disp( '该病例没有此牙齿 ');
 
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
        disp( '该病例没有此牙齿 ');
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
  