function [V,E,H] = png2poly(filename,smooth_iters,max_points)
  % PNG2POLY Read a png file and use it's alpha mask to define a polygon
  %
  % [V,E,H] = png2poly(filename)
  % [V,E,H] = png2poly(filename,smooth_iters,max_points)
  %
  % Inputs:
  %   filename  path to .png file
  %   smooth_iters  optional number of smoothing iterations to perform on output
  %     polygon, {0}
  %   max_points  optional number of maximum points along boundary, {inf}
  % Outputs
  %   V  #V by 2 list of polygon vertices, note that y-coordinates will be
  %     "flipped" with respect to the image if you think of the image rows as
  %     counting from 0 at the top left corner to size(im,1) at the bottom
  %     right corner
  %   E  #E by 2 list of polygon edge indices
  %   H  #H by 2 list of hole positions
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %   
  % See also: mask2poly, png2objandtga, png2mesh
  %

  
  % read in .png file with alpha layer
  [img,map,alpha] = imread(filename);
  % let shape be area with alpha > 128
  alpha = 255*( alpha > 128);
  % get polygon from alpha mask
  poly = mask2poly(alpha);

  if exist('smooth_iters','var')
    % handle smooth of each component separately, *before* finding hole positions
    for component_index = 1:size(poly,2)
      % flip y-coordinates because matlab displays images with reversed
      % y-direction, but displays meshes with normal y-direction
      poly(component_index).y = size(alpha,1)-poly(component_index).y;
      component_size = size(poly(component_index).x,1);
        % this will probably not just work for non closed loops
        i = [ ...
          1:component_size ...
          [component_size 1:(component_size-1)] ...
          1:component_size ];
        j = [
          [component_size 1:(component_size-1)] ...
          1:component_size ...
          1:component_size ];
        v = [ ...
          0.5*ones(1,component_size) ...
          0.5*ones(1,component_size) ...
          -1*ones(1,component_size)];
        % normalized laplacian matrix
        L = sparse(i,j,v,component_size,component_size);
        for iteration = 1:smooth_iters
          % simple 1.0 parameter
          poly(component_index).x = ...
            [1.0*L * poly(component_index).x + poly(component_index).x];
          poly(component_index).y = ...
            [1.0*L * poly(component_index).y + poly(component_index).y];
        end
    end
  end


  % total number of vertices
  n = size(cat(1,poly.x),1);
  while(n > max_points)
    % downsample each component
    for component_index = 1:size(poly,2)
      if(n < max_points)
        break;
      end
      if(exist('max_points'))
        % reduce number of points by half
        poly(component_index).x = poly(component_index).x(1:2:end);
        poly(component_index).y = poly(component_index).y(1:2:end);
        component_size = size(poly(component_index).x,1);
      end
      % compute new total number of points
      n = size(cat(1,poly.x),1);
    end
  end

  % convert poly struct to V,E,H 
  [V,E,H] = poly2VEH(poly);
end
