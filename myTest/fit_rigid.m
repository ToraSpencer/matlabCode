function [R,t] = fit_rigid(vers1,vers2)
  S = bsxfun(@minus,vers2,mean(vers2))'*bsxfun(@minus,vers1,mean(vers1));
  % Find closest rotation
  R = fit_rotation(S');

  % translation is just difference of center of mass (assuming uniform mass
  % distribution)
  t = mean(vers1)-mean(vers2)*R;

end
