function L = cotmatrix_embedded(vers, tris)
  % ��������laplace����
 
  % Inputs:
  %   V  #V x dim matrix of vertex coordinates
  %   F  #F x 3  matrix of indices of triangle corners
  % Outputs:
  %   L  sparse #V x #V matrix of cot weights 
 
   % ����Ƭ�ı߳�����
   edgeLens = [ ...
     sqrt(sum((vers(tris(:,2),:)-vers(tris(:,3),:)).^2,2)) ...      % ��1������ĶԱ�
     sqrt(sum((vers(tris(:,3),:)-vers(tris(:,1),:)).^2,2)) ...      % ��2������ĶԱ�
     sqrt(sum((vers(tris(:,1),:)-vers(tris(:,2),:)).^2,2)) ...
     ];
 
   L = cotmatrix_intrinsic(edgeLens, tris, size(vers,1));
end

