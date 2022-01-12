%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                   triangulate_polygon 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Divide a polygon into triangles using triangulation Matlab's built-in function.
% This works for any polygon (convex or non-convex).
%
% Usage
% =====
% triangles = triangulate_polygon(domainMesh,element_id)
%
% Input
% =====
% domainMesh : structure containing the polygonal mesh information
% element_id : element number
%
% Output
% ======
% connect    : connectivity of the triangles
%
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
% Jan. 31, 2020: triangulate using Matlab's triangulation built-in function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function connect = triangulate_polygon(domainMesh,element_id)

%   fprintf('******* Triangulating element %d\n',element_id);
  
  elem_nodes=domainMesh.connect{element_id,1};
%   element_id
%   elem_nodes
  elem_coords=domainMesh.coords(elem_nodes,:);
%   T = triangulation(polyshape(elem_coords));
  T = triangulation(polyshape(elem_coords,'KeepCollinearPoints',true,'Simplify',false));  
  num_triangles = size(T.ConnectivityList,1);
  connect = zeros(num_triangles,3);
  for e = 1:num_triangles
    connect(e,:) = elem_nodes(T.ConnectivityList(e,:));
  end
  
%   connect
  
%   figure
%   triplot(T);
 
end

