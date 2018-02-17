%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                      quad4_mesh 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% 4-node quadrilateral element mesh generator for rectangular domains with 
% option to create uniform or nonuniform meshes.
%
% Usage
% =====
% [coords,connect] = quad4_mesh(mesh_type,xmin,ymin,xmax,ymax,ndiv_x,ndiv_y)
%
% Inputs
% ======                               
% mesh_type   : 'uniform' or 'nonuniform'
%                                               
% xmin        : x-coordinate of the lower left of the rectangular domain
% ymin        : y-coordinate of the lower left of the rectangular domain
% xmax        : x-coordinate of the upper right of the rectangular domain
% ymax        : y-coordinate of the upper right of the rectangular domain
% ndiv_x      : number of sub-divisions in the x-direction
% ndiv_y      : number of sub-divisions in the y-direction
%
% Outputs
% =======        
% coords   : coordinates of the nodes [x1 y1; x2 y2; ... ; xn yn]                     
% connect  : element connectivity (array of size (nelems)x(4))
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [coords,connect] = ...
                    quad4mesh(mesh_type,xmin,ymin,xmax,ymax,ndiv_x,ndiv_y)
       
  % compute nodes in x- and y-directions and total number of nodes 
  nodes_x=ndiv_x+1; 
  nodes_y=ndiv_y+1; 
  numnd=nodes_x*nodes_y; % numnd is the number of nodes in the mesh
  nelems=ndiv_x*ndiv_y;

  % nodes are equi-spaced in x-direction; hx = x-spacing
  x=linspace(xmin,xmax,nodes_x);
  % nodes are equi-spaced in y-direction; hy = y-spacing
  y=linspace(ymin,ymax,nodes_y);

  coords=zeros(numnd,2);
  for i=1:nodes_x
    indices=i:nodes_x:(nodes_x*(nodes_y-1)+i);
    coords(indices,1)=x(i).*ones(nodes_y,1);
    coords(indices,2)=y(1:nodes_y)';
  end

  % connectivities
  connect=zeros(nelems,4);
  for j=1:ndiv_y
    for i=1:ndiv_x
      e=ndiv_x*(j-1)+i;
      node=nodes_x*(j-1)+i;
      connect(e,:)=[node,node+1,node+nodes_x+1,node+nodes_x];
    end
  end
  
  % make unstructured mesh
  mesh_size=min((xmax-xmin)/ndiv_x,(ymax-ymin)/ndiv_y);
  if strcmp(mesh_type,'nonuniform')
    for i=1:numnd
      h=rand(1);
      hc=round(10*h);
      xp=coords(i,1:2);
      [boundary,~]=...
                 checkboundarypoint(xp,xmin,ymin,xmax,ymax,nodes_x,nodes_y);
      if (boundary == 0) % interior node
        if (rem(hc,2) == 0) % even
          coords(i,1)=coords(i,1)+h*rand(1)*mesh_size*0.4;
          coords(i,2)=coords(i,2)+h*rand(1)*mesh_size*0.4;
        else % odd
          coords(i,1)=coords(i,1)-h*rand(1)*mesh_size*0.4;
          coords(i,2)=coords(i,2)-h*rand(1)*mesh_size*0.4;
        end
      elseif (boundary == 1) % bootom boundary
        if (rem(hc,2) == 0) % even
          coords(i,1)=coords(i,1)+h*rand(1)*mesh_size*0.4;
        else % odd
          coords(i,1)=coords(i,1)-h*rand(1)*mesh_size*0.4;
        end     
      elseif (boundary == 2) % left boundary
        if (rem(hc,2) == 0) % even
          coords(i,2)=coords(i,2)+h*rand(1)*mesh_size*0.4;
        else % odd
          coords(i,2)=coords(i,2)-h*rand(1)*mesh_size*0.4;
        end  
      elseif (boundary == 3) % right boundary
        if (rem(hc,2) == 0) % even
          coords(i,2)=coords(i,2)+h*rand(1)*mesh_size*0.4;
        else % odd
          coords(i,2)=coords(i,2)-h*rand(1)*mesh_size*0.4;
        end 
      elseif (boundary == 4) % upper boundary
         if (rem(hc,2) == 0) % even
          coords(i,1)=coords(i,1)+h*rand(1)*mesh_size*0.4;
         else % odd
          coords(i,1)=coords(i,1)-h*rand(1)*mesh_size*0.4;
         end  
      end % when boundary=5 (corner node), there is no correction to coord
    end
  end
  
  % plot mesh
  quadrilaterals=cell(nelems,1);
  for e=1:nelems
    quadrilaterals{e,1}=connect(e,:)';
  end 
  figure; 
  maxNumVertices = max(cellfun(@numel,quadrilaterals));
  padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
  elements = cellfun(padFunc,quadrilaterals,'UniformOutput',false);
  elements = vertcat(elements{:});
  patch('Faces',elements,'Vertices',coords,'FaceColor','w'); 
  axis equal; axis off; 
  
end

function [boundary,cornernodenum] = ...
                checkboundarypoint(point,xmin,ymin,xmax,ymax,nodes_x,nodes_y)

  % check if the "point" belongs to the boundary of a square/rectangular
  % domain or to a corner of the domain.
  %
  % boundary = 0, if the node is not on the boundary
  % boundary = 1, if the node is on the bottom boundary
  % boundary = 2, if the node is on the left boundary
  % boundary = 3, if the node is on the right boundary
  % boundary = 4, if the node is on the upper boundary
  % boundary = 5, if the node is in one corner of the boundary
  % cornernodenum : is the global node number of a corner node

  boundary = 0;
  cornernodenum = 0;
  if ((point(1) == xmin) && (point(2) == ymin))
      boundary = 5;
      cornernodenum = 1;
  elseif ((point(1) == xmax) && (point(2) == ymin))
      boundary = 5;
      cornernodenum = nodes_x;    
  elseif ((point(1) == xmax) && (point(2) == ymax))
      boundary = 5;
      cornernodenum = nodes_x*nodes_y; 
  elseif ((point(1) == xmin) && (point(2) == ymax))
      boundary = 5;
      cornernodenum = nodes_x*(nodes_y-1)+1; 
  elseif (point(2) == ymin)
      boundary = 1;
  elseif (point(1) == xmin)
      boundary = 2;
  elseif (point(1) == xmax)
      boundary = 3;
  elseif (point(2) == ymax)
      boundary = 4;
  end
end


