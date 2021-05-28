clc;
clear all;
load('movedRootTooth1.mat');
load('point1.mat');
load('point2.mat');
load('movedRootTooth2.mat');


%%
for i = 1:length(point2)
   [minValue,r]=mindis(point1,point2(i,:),1);
   row(i) = r;
end

point1 = point1(row,:);

  max_iter = 100;
  min_R_norm = 1e-5;
  min_t_norm = 1e-5;

  max_samples = min([size(point1,1); size(point2,1)]);

  dim = size(point1,2);
  bbd = norm(max([point1; point2]) - min([point1;point2]));

  R = eye(dim);
  t = zeros(1,dim);
  tt = repmat(t, size(point2, 1), 1);
  BRt = point2*R + tt;
 
  KDTA = createns(point1,'nsmethod','kdtree');
  KDTB = createns(point2,'nsmethod','kdtree');

  iter = 1;
  
%%
  while true
    prev_R = R;
    prev_t = t;

    IA = randperm(size(point1,1));
    IA = IA(1:max_samples);
    IB = randperm(size(point1,1));
    IB = IB(1:max_samples);

    KAB = knnsearch(KDTA,BRt(IB,:),'K',1);
    
    KBA = knnsearch(KDTB,bsxfun(@minus,point1(IA,:),t)*R','K',1);

    vers1 = [point1(KAB,:);point1(IA,:)];
    vers2 = [point2(IB,:);point2(KBA,:)];
    center1 = mean(vers1);
    center2 = mean(vers2);
    centerMat1 = repmat(center1, size(vers1, 1),1);
    centerMat2 = repmat(center2, size(vers2, 1),1);
    S = (vers2 - centerMat2)'*(vers1 - centerMat1);
    
 
      [su,~,sv] = svd(S');     %奇异值分解
        R = sv*su';
      if( det(R) < 0 )
        su(:,end) = -su(:,end);
        R = sv*su';
      end
 
    
    t = mean(vers1)-mean(vers2)*R;
    tt = repmat(t, size(point2, 1), 1);
    BRt = point2*R + tt;

    if iter >1 && norm(R-prev_R)<min_R_norm && norm(t-prev_t)<min_t_norm*bbd
      break;
    end
    
    iter = iter + 1;
   
    if iter > max_iter
      break;
    end
  end

%%
t = repmat(t, size(movedRootTooth1, 1), 1);
movedRootTooth22 = movedRootTooth1 * R + t;

OBJwriteVertices('第二次旋转平移后的带根标准牙.obj', movedRootTooth22);
disp('finished.');
