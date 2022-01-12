%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                    xVEMLab
%                          Source code  : not released
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                   quad4_mesh_with_crack
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% 4-node quadrilateral element mesh generator for rectangular domains with 
% option to create uniform or nonuniform meshes.
%
% Usage
% =====
% [coords,connect,enriched_nodes] = quad4_mesh(mesh_type,xmin,ymin,xmax,ymax,ndiv_x,ndiv_y,crack_tips_coords)
%
% Input
% =====                               
% mesh_type   : 'uniform' or 'nonuniform'
%                                               
% xmin        : x-coordinate of the lower left of the rectangular domain
% ymin        : y-coordinate of the lower left of the rectangular domain
% xmax        : x-coordinate of the upper right of the rectangular domain
% ymax        : y-coordinate of the upper right of the rectangular domain
% ndiv_x      : number of sub-divisions in the x-direction
% ndiv_y      : number of sub-divisions in the y-direction
% crack_tips_coords : coordinates of the crack tips
% pltoptions  : plotting options
%
% Output
% ======        
% coords   : coordinates of the nodes [x1 y1; x2 y2; ... ; xn yn]                     
% connect  : element connectivity (array of size (nelems)x(4))
% enriched_nodes : [...; node_number enrichment_type; ...]; (enrichment_type: 1 (Heaviside) or 2 (crack-tip))
% 
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Dec. 26, 2019: initial release (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coords,connect,enriched_nodes] = ...
                    quad4mesh_with_crack(mesh_type,xmin,ymin,xmax,ymax,ndiv_x,ndiv_y,...
                    crack_tips_coords,pltoptions)
       
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
  
  % node numbering
  if pltoptions.plot_node_numbers
    offset_x = 0.06*mesh_size; offset_y = 0.08*mesh_size;
    for i = 1:numnd
      idstring = sprintf('%d',i);
      text(coords(i,1)+offset_x, coords(i,2)+offset_y, idstring,'color',[1,0,0],'fontsize',pltoptions.fsize);
    end
  end
  
  % element numbering
  if pltoptions.plot_element_numbers
    for e = 1:nelems
      xe = coords(connect(e,:),1);
      ye = coords(connect(e,:),2);
      polyin = polyshape(xe,ye);
      [xc,yc] = centroid(polyin);
      xc = xc+pltoptions.element_number_offset_x;
      yc = yc+pltoptions.element_number_offset_y;
      idstring = sprintf('%d',e);
      text(xc, yc, idstring, 'color',[0,0,0],'fontsize',pltoptions.fsize+2);
    end 
  end
  
  % crack plot
  if pltoptions.plot_crack
    xcrack = crack_tips_coords(:,1);
    ycrack = crack_tips_coords(:,2);
    line(xcrack,ycrack,'Color','blue','LineStyle','-','LineWidth',pltoptions.linesize);
  end
  
  % ask for enriched nodes
  enriched_nodes=input('Enter enriched nodes in the format [...; node_number enrichment_type; ...] (enrichment_type: 1 (Heaviside) or 2 (crack-tip))\n');  
  
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


