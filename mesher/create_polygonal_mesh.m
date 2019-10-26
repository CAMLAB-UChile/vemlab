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
  %%%%%%%%%%%%%%%%%%           USER AREA         %%%%%%%%%%%%%%%%%%%%%%%%%
  
  % NOTE: 
  % To change the geometry of the customized domains, see instructions in the 
  % following functions that are below in this file:
  %    - RectangularDomain
  %    - WrenchDomain
  %    - PlateWithHoleDomain
  %
  
  N_elems=4000;
  mesh_filename='big_wrench_4000poly_elems.txt'; 
  
  % WARNING: DON'T CHANGE THE FOLLOWING TWO LINES IF YOU DON'T KNOW WHAT YOU ARE DOING!
  config=config_vemlab_mesher(opsystem,vemlab_root_dir,mesh_filename); % configure mesher
  mesh_file=[config.mesh_folder_location,mesh_filename];
  
  % USE ONLY ONE OF THE FOLLOWING PolyMesher FUNCTIONS FOR CUSTOMIZED DOMAINS
  
  % RECTANGULAR DOMAIN
%   [~,~,~,~,~]=PolyMesher(@RectangularDomain,N_elems,100,mesh_file);
  
  % WRENCH DOMAIN  
%   [~,~,~,~,~]=PolyMesher(@WrenchDomain,N_elems,100,mesh_file);  

  % BIG WRENCH DOMAIN  
  [~,~,~,~,~]=PolyMesher(@BigWrenchDomain,N_elems,100,mesh_file);  
  
  % PLATE WITH HOLE DOMAIN  
%   [~,~,~,~,~]=PolyMesher(@PlateWithHoleDomain,N_elems,100,mesh_file);  

  
  %%%%%%%%%%%%%%%%%%         END USER AREA         %%%%%%%%%%%%%%%%%%%%%%%  
  
end

%----------------- PolyMesher's INTERFACE FUNCTIONS ----------------------%
% NOTE:                                                                   
%   - These are customized functions prepared by A. Ortiz-Bernardin using   
%   the examples provided in PolyMesher's original source code.           
%   Only rectangular and wrench domains are implemented. Other domain     
%   implementations must be added below. 
%   - Adding new customized domains also requires updating the following 
%     functions:
%       read_mesh.m (located in folder vemlab/preprocessing)
%       PolyMesher2Veamy.m (located in folder mesher/PolyMesher)
% 
% Function's updates history
% ==========================
% Dec. 26, 2017: initial release
% May 15, 2018: add wrench domain     
% May 16, 2018: add plate with a hole domain                   
%-------------------------------------------------------------------------%

%% 
% RECTANGULAR DOMAIN
%
function [x] = RectangularDomain(Demand,Arg)
  BdBox = [0 1 0 1]; % [xmin xmax ymin ymax]: in this case, BdBox defines
                     % the dimensions of the rectangular domain. Update the 
                     % BdBox to change the dimensions.
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
    case('Boundary'); x = BoundaryRectangularDomain(Arg{:},BdBox); % added by A.O-B
    case('DomainType'); x = 'RectangularDomain'; % added by A.O-B
  end
  %----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
  function Dist = DistFnc(P,BdBox)
   Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));   
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
end

%% 
% WRENCH DOMAIN
%
function [x] = WrenchDomain(Demand,Arg)
  BdBox = [-0.3 2.5 -0.5 0.5]; % [xmin xmax ymin ymax]: in this case, BdBox
                               % defines a box within which the wrench domain
                               % will reside. The actual wrench domain is
                               % defined in function DistFnc below. 
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
    case('Boundary'); x = BoundaryWrenchDomain(Arg{:},BdBox); % added by A.O-B 
    case('DomainType'); x = 'WrenchDomain'; % added by A.O-B      
  end
  %----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
  function Dist = DistFnc(P,BdBox) 
    % Update this function to change the dimensions of the wrench domain;
    % if updated also update BdBox above.
    d1 = dLine(P,0,0.3,0,-0.3);
    d2 = dLine(P,0,-0.3,2,-0.5);
    d3 = dLine(P,2,-0.5,2,0.5);
    d4 = dLine(P,2,0.5,0,0.3);
    d5 = dCircle(P,0,0,0.3);
    d6 = dCircle(P,2,0,0.5);
    douter = dUnion(d6,dUnion(d5,...
             dIntersect(d4,dIntersect(d3,dIntersect(d2,d1)))));
    d7 = dCircle(P,0,0,0.175);
    d8 = dCircle(P,2,0,0.3);
    din = dUnion(d8,d7);
    Dist = dDiff(douter,din);
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
  function [BoundaryNodes] = BoundaryWrenchDomain(Node,BdBox)
    eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
    BoundaryNodes.RightCircle = ...
        find(abs(sqrt((Node(:,1)-2).^2+ Node(:,2).^2)-0.3)<eps);   
    BoundaryNodes.LeftHalfCircle = ...
        find(abs(max(sqrt(Node(:,1).^2+Node(:,2).^2)-0.175,Node(:,2)))<eps);    
  end 
end

%% 
% BIG WRENCH DOMAIN
%
function [x] = BigWrenchDomain(Demand,Arg)
  BdBox = [-8 63 -13 13]; % [xmin xmax ymin ymax]: in this case, BdBox
                               % defines a box within which the wrench domain
                               % will reside. The actual wrench domain is
                               % defined in function DistFnc below.    
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
    case('Boundary'); x = BoundaryWrenchDomain(Arg{:},BdBox); % added by A.O-B 
    case('DomainType'); x = 'WrenchDomain'; % added by A.O-B      
  end
  %----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
  function Dist = DistFnc(P,BdBox) 
    % Update this function to change the dimensions of the wrench domain;
    % if updated also update BdBox above.
    d1 = dLine(P,0,8,0,-8);
    d2 = dLine(P,0,-8,50,-13);
    d3 = dLine(P,50,-13,50,13);
    d4 = dLine(P,50,13,0,8);
    d5 = dCircle(P,0,0,8);
    d6 = dCircle(P,50,0,13);
    douter = dUnion(d6,dUnion(d5,...
             dIntersect(d4,dIntersect(d3,dIntersect(d2,d1)))));
    d7 = dCircle(P,0,0,4);
    d8 = dCircle(P,50,0,7);
    din = dUnion(d8,d7);
    Dist = dDiff(douter,din);
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
  function [BoundaryNodes] = BoundaryWrenchDomain(Node,BdBox)
    eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
    BoundaryNodes.RightCircle = ...
        find(abs(sqrt((Node(:,1)-50).^2+ Node(:,2).^2)-7)<eps);   
    BoundaryNodes.LeftHalfCircle = ...
        find(abs(max(sqrt(Node(:,1).^2+Node(:,2).^2)-4,Node(:,2)))<eps);    
  end 
end

%%
% PLATE WITH A HOLE DOMAIN
%
function [x] = PlateWithHoleDomain(Demand,Arg)
  BdBox = [0 5 -2 2]; % [xmin xmax ymin ymax]: in this case, BdBox defines
                      % the dimensions of the rectangular domain for the plate
                      % but WITHOUT the hole. The hole is added in the 
                      % DistFnc function below. Update the BdBox to change the 
                      % dimensions of the main plate domain.
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
    case('Boundary'); x = BoundaryPlateWithHoleDomain(Arg{:},BdBox); % added by A.O-B 
    case('DomainType'); x = 'PlateWithHoleDomain'; % added by A.O-B        
  end
 
  %----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
  function Dist = DistFnc(P,BdBox)
    % Update this function to change the dimensions of the plate 
    % with a hole domain    
    d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
    d2 = dCircle(P,0,0,BdBox(4)/2); % update this for the hole diameter (last argument)
    Dist = dDiff(d1,d2);
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
  function [BoundaryNodes] = BoundaryPlateWithHoleDomain(Node,BdBox)
    eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
    %BoundaryNodes.LeftCircle = find(abs(sqrt(Node(:,1).^2+Node(:,2).^2)-1.0)<eps);
    BoundaryNodes.Left = find(abs(Node(:,1)-BdBox(1))<eps);    
    BoundaryNodes.Right = find(abs(Node(:,1)-BdBox(2))<eps);
    MidRightFace = sqrt((Node(:,1)-BdBox(2)).^2+...
                        (Node(:,2)-(BdBox(3)+BdBox(4))/2).^2);
    [foo,MidRightFace] = sort(MidRightFace);  
    BoundaryNodes.MidNodeRightFace=MidRightFace(1);
  end  
end


