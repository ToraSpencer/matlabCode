clc;
clear all;
load('movedRootTooth1.mat');
load('point1.mat');
load('point2.mat');
load('movedRootTooth2.mat');


%% 1
for i = 1:length(point2)
   [minValue,r] = mindis(point1,point2(i,:),1);
   row(i) = r;
end


  point1 = point1(row,:);

  maxLoopCount = 100;
 

  sampleSize = min([size(point1,1); size(point2,1)]);
  
  KDTA = createns(point1,'nsmethod','kdtree');
  KDTB = createns(point2,'nsmethod','kdtree');

  loopCount = 1;
  
    % for output
    OBJwriteVertices('point1.obj', point1);
    OBJwriteVertices('point2.obj', point2);
  
    
  tempVers = [point1; point2];
  temp1 = max(tempVers);
  temp2 = min(tempVers);
 temp = temp1 - temp2;
  bbd = norm(temp);
 
  
  R = eye(3);
  t = zeros(1,3);
  tt = repmat(t, size(point2, 1), 1);
  BRt = point2 * R + tt;
 
%% 2
  while true
    prev_R = R;
    prev_t = t;

    IA = randperm(size(point1,1));
    IA = IA(1:sampleSize);
    IB = randperm(size(point1,1));
    IB = IB(1:sampleSize);

 


    tempb = BRt(IB,:);                       % 397个顶点坐标   
    tempa = point1(IA,:);
    tempa = tempa - repmat(t, length(point1),1);
    tempa = tempa*R';
    KAB = knnsearch(KDTA,tempb,'K',1);       % 列向量， 
    KBA = knnsearch(KDTB, tempa, 'K', 1);
   
    point1ab = point1(KAB,:);
    point1a = point1(IA,:);
    point2b = point2(IB,:);
    point2ba = point2(KBA,:);
    
    vers1 = [point1ab; point1a];
    vers2 = [point2b; point2ba];
    center1 = mean(vers1);
    center2 = mean(vers2);
    centerMat1 = repmat(center1, size(vers1, 1),1);
    centerMat2 = repmat(center2, size(vers2, 1),1);
    temp1 = (vers2 - centerMat2)';
    temp2 = (vers1 - centerMat1);
    S = temp1*temp2;
    
  [su,~,sv] = svd(S');     %奇异值分解
   R = sv*su';
  if( det(R) < 0 )
    su(:,end) = -su(:,end);
    R = sv*su';
  end
 
    
    t = mean(vers1)-mean(vers2)*R;
    tt = repmat(t, size(point2, 1), 1);
    BRt = point2*R + tt;

    if loopCount >1 && norm(R-prev_R)< 1e-5 && norm(t-prev_t) < 1e-5 * bbd
      break;
    end
    
    loopCount = loopCount + 1;
   
    if loopCount > maxLoopCount
      break;
    end
  end

%% 3
t = repmat(t, size(movedRootTooth1, 1), 1);
movedRootTooth22 = movedRootTooth1 * R + t;

OBJwriteVertices('第二次旋转平移后的带根标准牙.obj', movedRootTooth22);
disp('finished.');
