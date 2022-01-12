%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                      gunelvemesh 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Guñelve style mesh generator for rectangular domains.
%
% Usage
% =====
% [coords,connect,hx,hy] = gunelvemesh(mesh_type,xmin,ymin,xmax,ymax,ndiv_x,ndiv_y)
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
%
% Output
% ======        
% coords   : coordinates of the nodes [x1 y1; x2 y2; ... ; xn yn]                     
% connect  : element connectivity (array of size (nelems)x(4))
% hx       : size of an element in x-direction (rectangular domain!)
% hy       : size of an element in y-direction (rectangular domain!)
% 
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Apr. 26, 2021: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coords,connect,hx,hy] = ...
                    gunelvemesh(xmin,ymin,xmax,ymax,npattern_x,npattern_y)
                  
 
 close all;
 clc;
 
 A=xmax-xmin;
 B=ymax-ymin;
 a=A/npattern_x;
 b=B/npattern_y;
 hx=a/2; % size of an element in x-direction (rectangular domain!)
 hy=b/2; % size of an element in y-direction (rectangular domain!)
 
 % nodal coordinates
 if npattern_y==1
   numnd=2*npattern_x+1+3*npattern_x+1+2*npattern_x+1+14*npattern_x;
 elseif npattern_y>1
   numnd=(2*npattern_x+1+3*npattern_x+1+2*npattern_x+1+14*npattern_x)*npattern_y-(2*npattern_x+1)*(npattern_y-1);
 else
   throw_error('Error in gunelvemesh.m: npattern_y');
 end
 coords=zeros(numnd,2);
 
 % connectivity
 nelems=npattern_x*npattern_y*5;
 connect=cell(nelems,1); 
 
 
 %%%%%
 % form the first row of the mesh
 %%%%  
 
 % First pattern = master pattern
 
 % Master pattern located on the left bottom corner.
 % Nodes of the boundary of the master element: counterclock-wise starting
 % from the left bottom corner
 x1=xmin; y1=ymin;
 x2=x1+a/2; y2=y1;
 x3=x1+a; y3=y1;
 x4=x1+a; y4=y1+b/2;
 x5=x1+a; y5=y1+b;
 x6=x1+a/2; y6=y1+b;
 x7=x1; y7=y1+b;
 x8=x1; y8=y1+b/2;
 % nodes of the eight-pointed star
 %x9=x1+a/8; y9=y1+3*b/8;
 x9=x1+2*a/8; y9=y1+3*b/8; 
 x10=x1+3*a/8; y10=y1+3*b/8;
 %x11=x1+3*a/8; y11=y1+b/8;
 x11=x1+3*a/8; y11=y1+2*b/8; 
 %x12=x1+a/2; y12=y1+b/4;
 x12=x1+a/2; y12=y1+3*b/8; 
 %x13=x1+5*a/8; y13=y1+b/8;
 x13=x1+5*a/8; y13=y1+2*b/8; 
 x14=x1+5*a/8; y14=y1+3*b/8;
 %x15=x1+7*a/8; y15=y1+3*b/8;
 x15=x1+6*a/8; y15=y1+3*b/8; 
 %x16=x1+6*a/8; y16=y1+b/2;
 x16=x1+5*a/8; y16=y1+b/2; 
 %x17=x1+7*a/8; y17=y1+5*b/8;
 x17=x1+6*a/8; y17=y1+5*b/8;
 x18=x1+5*a/8; y18=y1+5*b/8;
 %x19=x1+5*a/8; y19=y1+7*b/8;
 x19=x1+5*a/8; y19=y1+6*b/8; 
 %x20=x1+a/2; y20=y1+6*b/8;
 x20=x1+a/2; y20=y1+5*b/8; 
 %x21=x1+3*a/8; y21=y1+7*b/8;
 x21=x1+3*a/8; y21=y1+6*b/8; 
 x22=x1+3*a/8; y22=y1+5*b/8;
 %x23=x1+a/8; y23=y1+5*b/8;
 x23=x1+2*a/8; y23=y1+5*b/8;
 %x24=x1+a/4; y24=y1+b/2; 
 x24=x1+3*a/8; y24=y1+b/2;
 
 % nodal coordinates of the master pattern
 x_master_pattern=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24]';
 y_master_pattern=[y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21,y22,y23,y24]';
 coords(1:24,1)=x_master_pattern;
 coords(1:24,2)=y_master_pattern; 
 nodes_master_pattern=1:24;
 nodes_master_to_copy_x=[2,3,4,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
 nodes_master_to_copy_y=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]; 
 nodes_master_to_copy_xy=[4,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
 
 % LINEAR ZOOM EFFECT: scales the master pattern (COMMENT IF NOT REQUIRED)
 %    Scales from a square/rectangular domain to another square/rectangular
 %    domain having both the same center
 % TRANSFORMATION: T(x,y)=(a1+a2*(x-xc),a3+a4*(y-yc)
 %     where (xc,yc) is the center of the square/rectangular domain
 %
 xx1=x18;
 yy1=y18;
 xx2=xx1+0.05*a;
 yy2=yy1+0.05*b;
 xx3=x10;
 yy3=y10;
 xx4=xx3-0.05*a;
 yy4=yy3-0.05*b;
 xc=xmin+A/2;
 yc=ymin+B/2;
 a1=xx2-(xx2-xx4)*(xx1-xc)/(xx1-xx3);
 a2=(xx2-xx4)/(xx1-xx3);
 a3=yy2-(yy2-yy4)*(yy1-yc)/(yy1-yy3);
 a4=(yy2-yy4)/(yy1-yy3);
 coords(9:end,:)=[a1+a2*(coords(9:end,1)-xc), a3+a4*(coords(9:end,2)-yc)];
 
 % master pattern connectivity
 connect{1,1}=[1,2,12,11,10,9,24,8];
 connect{2,1}=[2,3,4,16,15,14,13,12]; 
 connect{3,1}=[16,4,5,6,20,19,18,17]; 
 connect{4,1}=[8,24,23,22,21,20,6,7];  
 connect{5,1}=[9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];  
 pattern_nelems=5; 
 
 row{1,1}=nodes_master_pattern; % nodes of the first pattern
 
 nodes_copied_pattern=zeros(1,length(nodes_master_pattern));

 % complete first row 
 for j=2:npattern_x
   nodes_to_copy=row{1,j-1}(nodes_master_to_copy_x);
   size1=length(nodes_to_copy);
   range1=(nodes_to_copy(size1)+1):(nodes_to_copy(size1)+21);
   nodes_copied_pattern(1)=row{1,j-1}(nodes_master_pattern(3));
   nodes_copied_pattern(7)=row{1,j-1}(nodes_master_pattern(5)); 
   nodes_copied_pattern(8)=row{1,j-1}(nodes_master_pattern(4));  
   nodes_copied_pattern(nodes_master_to_copy_x)=range1; 
   coords(range1,1)=coords(nodes_to_copy,1)+a;
   coords(range1,2)=coords(nodes_to_copy,2); 
   connect{5*j-4,1}=nodes_copied_pattern(connect{1,1});
   connect{5*j-3,1}=nodes_copied_pattern(connect{2,1}); 
   connect{5*j-2,1}=nodes_copied_pattern(connect{3,1}); 
   connect{5*j-1,1}=nodes_copied_pattern(connect{4,1});  
   connect{5*j,1}=nodes_copied_pattern(connect{5,1});  
   row{1,j}=nodes_copied_pattern; % nodes of the j-th pattern of the first row   
 end

 % following rows
 for i=2:npattern_y
   last_element_previous_row=5*npattern_x*(i-1);
   for j=1:npattern_x
     if j==1
       nodes_to_copy=row{i-1,1}(nodes_master_to_copy_y);
       last_node_previous_row=row{i-1,npattern_x}(end);
       range1=(last_node_previous_row+1):(last_node_previous_row+21);
       nodes_copied_pattern(1)=row{i-1,1}(nodes_master_pattern(7));
       nodes_copied_pattern(2)=row{i-1,1}(nodes_master_pattern(6)); 
       nodes_copied_pattern(3)=row{i-1,1}(nodes_master_pattern(5)); 
       nodes_copied_pattern(nodes_master_to_copy_y)=range1; 
       coords(range1,1)=coords(nodes_to_copy,1);  
       coords(range1,2)=coords(nodes_to_copy,2)+b;
       connect{last_element_previous_row+5*j-4,1}=nodes_copied_pattern(connect{1,1});
       connect{last_element_previous_row+5*j-3,1}=nodes_copied_pattern(connect{2,1}); 
       connect{last_element_previous_row+5*j-2,1}=nodes_copied_pattern(connect{3,1}); 
       connect{last_element_previous_row+5*j-1,1}=nodes_copied_pattern(connect{4,1});  
       connect{last_element_previous_row+5*j,1}=nodes_copied_pattern(connect{5,1});   
       row{i,1}=nodes_copied_pattern; % nodes of the first pattern of the i-th row          
     else
       nodes_to_copy=row{i,j-1}(nodes_master_to_copy_xy); 
       size1=length(nodes_to_copy);
       range1=(nodes_to_copy(size1)+1):(nodes_to_copy(size1)+19);
       nodes_copied_pattern(1)=row{i,j-1}(nodes_master_pattern(3));
       nodes_copied_pattern(2)=row{i-1,j}(nodes_master_pattern(6)); 
       nodes_copied_pattern(3)=row{i-1,j}(nodes_master_pattern(5));  
       nodes_copied_pattern(7)=row{i,j-1}(nodes_master_pattern(5));
       nodes_copied_pattern(8)=row{i,j-1}(nodes_master_pattern(4)); 
       nodes_copied_pattern(nodes_master_to_copy_xy)=range1;
       coords(range1,1)=coords(nodes_to_copy,1)+a;
       coords(range1,2)=coords(nodes_to_copy,2);  
       connect{last_element_previous_row+5*j-4,1}=nodes_copied_pattern(connect{1,1});
       connect{last_element_previous_row+5*j-3,1}=nodes_copied_pattern(connect{2,1}); 
       connect{last_element_previous_row+5*j-2,1}=nodes_copied_pattern(connect{3,1}); 
       connect{last_element_previous_row+5*j-1,1}=nodes_copied_pattern(connect{4,1});  
       connect{last_element_previous_row+5*j,1}=nodes_copied_pattern(connect{5,1});  
       row{i,j}=nodes_copied_pattern; % nodes of the j-th pattern of the i-th row    
     end
   end
 end
 
%   % plot mesh
%   polygons=cell(nelems,1);
%   for e=1:nelems
%     polygons{e,1}=connect{e,1}';
%   end 
%   figure; 
%   maxNumVertices = max(cellfun(@numel,polygons));
%   padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
%   elements = cellfun(padFunc,polygons,'UniformOutput',false);
%   elements = vertcat(elements{:});
%   patch('Faces',elements,'Vertices',coords,'FaceColor','w'); 
%   axis equal; axis off; 
  
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


