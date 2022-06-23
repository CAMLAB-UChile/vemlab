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
  
  N_elems=2500;
  mesh_filename='plate_with_hole_2500poly_elems.txt';
  
  % WARNING: DON'T CHANGE THE FOLLOWING TWO LINES IF YOU DON'T KNOW WHAT YOU ARE DOING!
  config=config_vemlab_mesher(opsystem,vemlab_root_dir,mesh_filename); % configure mesher
  mesh_file=[config.mesh_folder_location,mesh_filename];
  
  % USE ONLY ONE OF THE FOLLOWING PolyMesher FUNCTIONS FOR CUSTOMIZED DOMAINS
  
  % RECTANGULAR DOMAIN
%   [Node,Element,NElem,BoundaryNodes,DomainType,BdBox,~,~,~]=...
%                      PolyMesher(@RectangularDomain,N_elems,1000,mesh_file);
  
  % WRENCH DOMAIN  
%   [Node,Element,NElem,BoundaryNodes,DomainType,BdBox,~,~,~]=...
%                      PolyMesher(@WrenchDomain,N_elems,100,mesh_file);  

  % BIG WRENCH DOMAIN  
%   [Node,Element,NElem,BoundaryNodes,DomainType,BdBox,~,~,~]=...
%                      PolyMesher(@BigWrenchDomain,N_elems,100,mesh_file);  
  
  % PLATE WITH HOLE DOMAIN  
  [Node,Element,NElem,BoundaryNodes,DomainType,BdBox,~,~,~]=...
                     PolyMesher(@PlateWithHoleDomain,N_elems,100,mesh_file);  

  % distort the mesh (comment if not required)
%   b=3.4;
%   a=0;
%   Node(:,1)=Node(:,1)+0.3*Node(:,1).*sin(b*Node(:,1)+a).*sin(b*Node(:,2)+a).*(BdBox(2)-Node(:,1))./BdBox(2);
%   Node(:,2)=Node(:,2)+0.3*Node(:,2).*sin(b*Node(:,2)+a).*sin(b*Node(:,1)+a).*(BdBox(4)-Node(:,2))./BdBox(4);

  %Plot mesh to a VEMLab mesh format
  PolyMesher2VEMLab(Node,Element,NElem,BoundaryNodes,mesh_file,DomainType,BdBox);      
 
  % plot mesh
  config.plot_mesh = 'yes'; % 'yes' or 'no'
  config.plot_node_numbers = 'no'; % 'yes' or 'no'
  config.plot_element_numbers = 'no'; % 'yes' or 'no'  
  config.plot_axis = 'no'; % 'yes' or 'no'
  %%%%%%%%%%%%%%%%%%         END USER AREA         %%%%%%%%%%%%%%%%%%%%%%%  
  plot_mesh(Element,Node,config);  
  
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
function [x,Arg] = RectangularDomain(Demand,Arg)
  BdBox = [0 2 0 2]; % [xmin xmax ymin ymax]: in this case, BdBox defines
                     % the dimensions of the rectangular domain. Update the 
                     % BdBox to change the dimensions.
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
    case('Boundary'); [x,Arg] = BoundaryRectangularDomain(Arg{:},BdBox); % added by A.O-B
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
  function [BoundaryNodes,Node] = BoundaryRectangularDomain(Node,BdBox)
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
    % fix boundary nodes coordinates to make them lie exactly on the BdBox
    Node(BoundaryNodes.bottom,2)=BdBox(3);    
    Node(BoundaryNodes.top,2)=BdBox(4); 
    Node(BoundaryNodes.left,1)=BdBox(1);      
    Node(BoundaryNodes.right,1)=BdBox(2);    
  end  
end

%% 
% WRENCH DOMAIN
%
function [x,Arg] = WrenchDomain(Demand,Arg)
  BdBox = [-0.3 2.5 -0.5 0.5]; % [xmin xmax ymin ymax]: in this case, BdBox
                               % defines a box within which the wrench domain
                               % will reside. The actual wrench domain is
                               % defined in function DistFnc below. 
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
    case('Boundary'); [x,Arg] = BoundaryWrenchDomain(Arg{:},BdBox); % added by A.O-B 
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
  function [BoundaryNodes,Node] = BoundaryWrenchDomain(Node,BdBox)
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
function [x,Arg] = BigWrenchDomain(Demand,Arg)
  BdBox = [-8 63 -13 13]; % [xmin xmax ymin ymax]: in this case, BdBox
                               % defines a box within which the wrench domain
                               % will reside. The actual wrench domain is
                               % defined in function DistFnc below.    
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
    case('Boundary'); [x,Arg] = BoundaryWrenchDomain(Arg{:},BdBox); % added by A.O-B 
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
  function [BoundaryNodes,Node] = BoundaryWrenchDomain(Node,BdBox)
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
function [x,Arg] = PlateWithHoleDomain(Demand,Arg)
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
    case('Boundary'); [x,Arg] = BoundaryPlateWithHoleDomain(Arg{:},BdBox); % added by A.O-B 
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
  function [BoundaryNodes,Node] = BoundaryPlateWithHoleDomain(Node,BdBox)
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

function plot_mesh(connect,coords,config)

  if strcmp(config.plot_mesh,'yes')
    nelems=length(connect);
    polygons=cell(nelems,1);
    for e=1:nelems
      polygons{e,1}=connect{e,1}';
    end 
    figure; 
    maxNumVertices = max(cellfun(@numel,polygons));
    padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
    elements = cellfun(padFunc,polygons,'UniformOutput',false);
    elements = vertcat(elements{:});
    patch('Faces',elements,'Vertices',coords,'FaceColor','w'); 
    axis equal; 
    if strcmp(config.plot_axis,'yes')
     axis on;
    else
      axis off;
    end

    hold on;

    % set the font size for the node and element numbering
    if nelems <= 50
      fsize = 10;
    else
      fsize = 6;
    end

    % node numbering
    numnod=size(coords,1);
    mesh_size=min_edge_size(coords,connect,nelems);
    if strcmp(config.plot_node_numbers,'yes')
      offset_x = 0.06*mesh_size; offset_y = 0.08*mesh_size;
      for i = 1:numnod
        idstring = sprintf('%d',i);
        text(coords(i,1)+offset_x, coords(i,2)+offset_y, idstring,'color',[1,0,0],'fontsize',fsize);
      end
    end

    % element numbering for QUAD4
    if strcmp(config.plot_element_numbers,'yes')
      for e = 1:length(connect(:,1))
        xc = mean(coords(connect{e},1)); % element number is plotted in the
        yc = mean(coords(connect{e},2)); % center of the rectangle
        idstring = sprintf('%d',e);
        text(xc, yc, idstring, 'color',[0,0,0],'fontsize',fsize);
      end
    end  
  end
  
end

function h_min = min_edge_size(Node,Element,NElem)
  h_min=1e1000;
  for e=1:NElem
    NVertex = length(Element{e});
    nodes = Element{e}(1:NVertex);
    verts=Node(nodes,:); 

    mysize=size(verts,1);
    mydist=zeros(mysize*mysize,1);

    % This procedure also compute the distance between a node with itself (0); however, it
    % doesn't matter since at the end the maximum distance is selected
    r = 1;
    for node_i = 1:mysize
      for node_j = 1:mysize
        % shifts coordinates with respect to "node i" in the connectivity
        shifted_coords= verts-repmat(verts(node_i,:),mysize,1);
        % compute the distances wrt to "node_i" in the element connectivity
        for node = 1:(mysize-1)
          mydist(r) = norm(shifted_coords(node,:));
          r = r + 1;
        end
      end  
    end

    % minimum distance from the first node to any node in the element connectivity 
    h_size=min(mydist);

    if h_size<h_min
      h_min=h_size;
    end   
  end
end



