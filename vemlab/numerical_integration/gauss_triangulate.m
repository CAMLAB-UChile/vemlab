%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                   gauss_triangulate
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Divide a polygon into triangles using the center of the polygon.
%
% Usage
% =====
% [gauss_points,gauss_weights] = ...
%     gauss_triangulate(elem_coords,elem_connect,gauss_order,plot_triangulation)
%
% Input
% =====
%   elem_coords : coordinates of 1 polygonal element (1x2 array)
%   elem_connect: connectivity of 1 polygonal element 
%                (1xN array, N number of nodes that form the polygonal element)
%   gauss_order : order for the quadrature for each of the triangles in the
%                 triangulation of the polygonal element (1, 2, 3, 4 or 5)
%   plot_triangulation : 'yes' or 'no' to plot the triangulation on the screen
%
% Output
% ======
% The following two arrays that contain the Gauss quadratures for the 
% polygonal element using triangulation:
%
%   gauss_points(num_gauss_points*num_triangles,2)
%   gauss_weights(num_gauss_points*num_triangles,1)
%
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gauss_points,gauss_weights] = ...
      gauss_triangulate(elem_coords,elem_connect,gauss_order,plot_triangulation)

  num_nodes=length(elem_connect);
  num_edges=num_nodes;

  %% Center of polygon
  center=[sum(elem_coords(:,1)),sum(elem_coords(:,2))]/num_nodes;
  coords=[elem_coords;elem_coords(1,:)];

  %% Triangulate using the center of the polygon, and generate Gauss points on each triangle of the triangulation
  num_gauss_points=num_gauss_points_T3(gauss_order);
  gauss_points=zeros(num_gauss_points*num_edges,2);
  gauss_weights=zeros(num_gauss_points*num_edges,1);
  if strcmp(plot_triangulation,'yes'), hold on; end
  for i=1:num_edges
    triangle=[coords(i,1),coords(i,2);coords(i+1,1),coords(i+1,2);center(1),center(2)];    
    [x,w]=gauss_points_T3(gauss_order,triangle);
    range=(num_gauss_points*i-num_gauss_points+1):(num_gauss_points*i);
    gauss_points(range,:)=x;
    gauss_weights(range)=w;
    if strcmp(plot_triangulation,'yes')
      for a=1:3
        for b=1:3
          plot([triangle(a,1);triangle(b,1)],[triangle(a,2);triangle(b,2)]);
        end
      end
    end
  end
  
  %% Plot Gauss points
  if strcmp(plot_triangulation,'yes')
      for i=1:length(gauss_points)
          plot(gauss_points(i,1),gauss_points(i,2),'x');
      end
  end

end
