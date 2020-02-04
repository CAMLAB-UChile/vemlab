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
% Divide a polygon into triangles using the polygon's vertices.
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function connect = triangulate_polygon_old(domainMesh,element_id)
  elem_nodes=domainMesh.connect(element_id,:);
  elem_num_nodes=length(elem_nodes{1});
  num_triangles=elem_num_nodes-2;
  connect=zeros(num_triangles,3);
  for elem_i=1:num_triangles
    connect(elem_i,:)=[elem_nodes{1}(1),elem_nodes{1}(elem_i+1),elem_nodes{1}(elem_i+2)];
  end  
end

