clc
close all;
clear all;

addpath('mesh process');
addpath('mixFE');  


%%
load('N0.mat');
load('N1.mat');
load('omega.mat');
load('mergedToothVers.mat');
load('mergeVers.mat');
load('newTris.mat');
load('bi_L.mat');
load('bi_U.mat');
load('bi_P.mat' );
load('bi_Q.mat' );
load('bi_R.mat' );
load('bi_S.mat' );
load('bi_M.mat');
load('patientTransform.mat');
load('index_tooth_change.mat');

  finalVers = mergedToothVers;
for i = 1:length(mergeVers)
   [minValue,r]=mindis(patientTransform,mergeVers(i,:),1);
   finalVers(index_tooth_change(i),:) = patientTransform(r,:);
end
  
  BZ1 = zeros(size(finalVers,1),3);
  all = [N0 omega]; 
  b0 = finalVers(N0,:);
  f1 = finalVers(N1,:);
  bn = zeros(size([N0 omega],2),3) ;
  
  for coord_index = 1:3
    rhs_Dx = -bi_S(all,N0)*b0(:,coord_index) - bi_S(all,N1)*f1(:,coord_index);
    rhs_Dy = zeros(size(omega,2),1);
    rhs = [ rhs_Dx; rhs_Dy];   
    sol = bi_Q*(bi_U\(bi_L\(bi_P*(bi_R\rhs))));
    ny = size(all,2);
    finalVers(omega,coord_index) = sol(ny+1:end);   
  end
  writeOBJ('×îÖÕ½á¹û.obj', finalVers, newTris);
  disp('finished');
  

