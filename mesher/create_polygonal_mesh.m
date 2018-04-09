function create_polygonal_mesh
  close all;
  clc;
  
  %% ADD ALL VEMLAB FOLDERS TO THE PATH
  opsystem=computer;
  is_Windows = strcmp(opsystem,'PCWIN') || strcmp(opsystem,'PCWIN64');
  is_Linux = strcmp(opsystem,'GLNX86') || strcmp(opsystem,'GLNXA64');   
  if is_Windows
    cd ..\; vemlab_root_dir=setpath;
  elseif is_Linux
    cd ../; vemlab_root_dir=setpath;
  end   

  %%
  %%%%%%%%%%%%%%%%%%           USER INPUT DATA         %%%%%%%%%%%%%%%%%%%%%%%%%
  
  % NOTE: you can change the rectangular domain dimensions by changing BdBox in
  % the RectangularDomain function below
  
  N_elems=2000;
  mesh_filename='square_plate_poisson2d_2000poly_elems.txt'; 
  
  %%%%%%%%%%%%%%%%%%%        END USER INPUT DATA       %%%%%%%%%%%%%%%%%%%%%%%%%   
  
  %%  
  config=config_vemlab_mesher(opsystem,vemlab_root_dir,mesh_filename); % configure mesher
  mesh_file=[config.mesh_folder_location,mesh_filename];
  [~,~,~,~,~]=PolyMesher(@RectangularDomain,N_elems,100,mesh_file);
  
end

%----------------- PolyMesher's INTERFACE FUNCTIONS ----------------------%
% NOTE:                                                                   %
%   These are customized functions prepared by A. Ortiz-Bernardin using   %
%   the examples provided in PolyMesher's original source code.           %
%   Only rectangular domain is implemented. Other domain implementations  %
%   must be added here.                                                   %
%                                                                         %
%   Dated: Dec. 26, 2017                                                  %
%-------------------------------------------------------------------------%
function [x] = RectangularDomain(Demand,Arg)
  BdBox = [0 1 0 1]; % [xmin xmax ymin ymax]
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
    case('Boundary'); x = BoundaryRectangularDomain(Arg{:},BdBox); % added by A.O-B
  end
end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  x1 = BdBox(1); x2 = BdBox(2); y1 = BdBox(3); y2 = BdBox(4);
  d = [x1-P(:,1), P(:,1)-x2, y1-P(:,2), P(:,2)-y2];
  d = [d,max(d,[],2)];  
  Dist = d;
end
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  x = {[],[]};
end
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
end
%----------------------------------------------------- GET BOUNDARY (added by A.O-B)
function [BoundaryNodes] = BoundaryRectangularDomain(Node,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  BoundaryNodes.bottom = find(abs(Node(:,2)-BdBox(3))<eps);   
  BoundaryNodes.top = find(abs(Node(:,2)-BdBox(4))<eps);  
  BoundaryNodes.left = find(abs(Node(:,1)-BdBox(1))<eps);
  BoundaryNodes.right = find(abs(Node(:,1)-BdBox(2))<eps);  
  BoundaryNodes.all = [BoundaryNodes.bottom;BoundaryNodes.top;...
                       BoundaryNodes.left;BoundaryNodes.right];
  BoundaryNodes.all=unique(BoundaryNodes.all); % delete repeated nodes from boundary_nodes.all ... possibly corner nodes
  blcorner_index=Node(BoundaryNodes.bottom,1)==min(Node(BoundaryNodes.bottom,1));
  BoundaryNodes.blcorner = BoundaryNodes.bottom(find(blcorner_index==1));
  trcorner_index=Node(BoundaryNodes.top,1)==max(Node(BoundaryNodes.top,1));
  BoundaryNodes.trcorner = BoundaryNodes.top(find(trcorner_index==1));
end
%-------------------------------------------------------------------------%