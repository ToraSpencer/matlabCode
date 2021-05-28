clc;
close all;
clear all;
 
%%
load('patientTooth.mat');
load('triIdxInCutPatient.mat');
load('raw_edges_list.mat');

tris = patientTooth.face(triIdxInCutPatient,:);

trisPlus = cat(2,tris,tris(:,1));  %   
R = repelem(trisPlus,1,[1 2 2 1]); %  ����Ƭ�е�[x,y,z]��дΪ[x,y,y,z,z,x];

edgeListTemp = cell2mat(cellfun(@(x) reshape(x,[2,3])',num2cell(R,2),'un',0));  % һ�б�����У�[x,y,y,z,z,x]�� [x,y; y,z; z,x]

edgeList = sort(edgeListTemp,2);   % ÿһ���������������������ǰС���
    
 disp('finished.');
 
 