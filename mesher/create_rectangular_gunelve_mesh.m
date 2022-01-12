function create_rectangular_gunelve_mesh
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
  
  xmin=0; xmax=2; 
  ymin=0; ymax=2;
  npattern_x=1; % number of patterns along x coordinate
  npattern_y=1; % number of patterns along y coordinate
  mesh_filename='colliding_flow_gunelvemesh_1x1_TEST.txt'; 
  
  %%%%%%%%%%%%%%%%%%%        END USER INPUT DATA       %%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%
  % configure mesher
  config=config_vemlab_mesher(opsystem,vemlab_root_dir,mesh_filename);
  
  % run quad4_mesh
  [coords,connect,hx,hy]=gunelvemesh(xmin,ymin,xmax,ymax,npattern_x,npattern_y);
  
  % get boundary nodes
  BdBox=[xmin,ymin;xmax,ymax];   
  BoundaryNodes=BoundaryRectangularDomain(coords,BdBox);
  
  % distmesh2VEMLab: write mesh to a VEMLab mesh format
  mesh_file=[config.mesh_folder_location,mesh_filename];  
  gunelvemesh2VEMLab(coords,connect,size(connect,1),BoundaryNodes,mesh_file);   
  
  % plot mesh
  config.plot_mesh = 'yes'; % 'yes' or 'no'
  config.plot_mesh_linewidth = 1.1; % config.plot_mesh_nodesize
  config.plot_node_numbers = 'no'; % 'yes' or 'no'
  config.plot_element_numbers = 'no'; % 'yes' or 'no'  
  config.plot_mesh_nodes = 'yes'; % 'yes' or 'no'
  config.plot_mesh_nodesize = 8;
  config.plot_axis = 'no'; % 'yes' or 'no'
  config.offset_hx = hx/2;
  config.offset_hy = hy/2;
  plot_gunelve_mesh(connect,coords,config);
 
end

%------------------- quad4mesh's INTERFACE FUNCTIONS ----------------------%
%   Only rectangular domain is implemented. Other domain implementations  %
%   must be added here.                                                   %
%                                                                         %
%   Dated: Dec. 31, 2017                                                  %
%-------------------------------------------------------------------------%

%---------------------------------------------------- GET BOUNDARY
function [BoundaryNodes] = BoundaryRectangularDomain(Node,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  BoundaryNodes.bottom = find(abs(Node(:,2)-BdBox(3))<eps);   
  BoundaryNodes.top = find(abs(Node(:,2)-BdBox(4))<eps);  
  BoundaryNodes.left = find(abs(Node(:,1)-BdBox(1))<eps);
  BoundaryNodes.right = find(abs(Node(:,1)-BdBox(2))<eps); 
  %-----------------------------------------------------------------------------
  % Nodes on each of the boundaries must be ordered in the physical direction,
  % either ascending or descending. This is required for applying the Neumann
  % boundary conditions --- on each segment of the Neumann boundary, we need 
  % the start and ending node of each segment. Another better approach would be 
  % to have available the connectivity of the boundaries, which we don't have. 
  % So, we order the boundaries using the coordinates along each of the four 
  % boundaries.
  %-----------------------------------------------------------------------------
  % order bottom boundary
  [~,ordered_indices]=sort(Node(BoundaryNodes.bottom,1),'ascend');
  BoundaryNodes.bottom=BoundaryNodes.bottom(ordered_indices); 
  % order top boundary
  [~,ordered_indices]=sort(Node(BoundaryNodes.top,1),'ascend');
  BoundaryNodes.top=BoundaryNodes.top(ordered_indices);   
  % order left boundary
  [~,ordered_indices]=sort(Node(BoundaryNodes.left,2),'ascend');
  BoundaryNodes.left=BoundaryNodes.left(ordered_indices);   
  % order right boundary
  [~,ordered_indices]=sort(Node(BoundaryNodes.right,2),'ascend');
  BoundaryNodes.right=BoundaryNodes.right(ordered_indices);  
  %-----------------------------------------------------------------------------
  BoundaryNodes.all = [BoundaryNodes.bottom;BoundaryNodes.top;...
                       BoundaryNodes.left;BoundaryNodes.right]; 
  % delete repeated nodes from boundary_nodes.all ... possibly corner nodes                   
  BoundaryNodes.all=unique(BoundaryNodes.all);   
  blcorner_index=Node(BoundaryNodes.bottom,1)==min(Node(BoundaryNodes.bottom,1));
  
  BoundaryNodes.blcorner = BoundaryNodes.bottom(find(blcorner_index==1));
  trcorner_index=Node(BoundaryNodes.top,1)==max(Node(BoundaryNodes.top,1));
  BoundaryNodes.trcorner = BoundaryNodes.top(find(trcorner_index==1));
end
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% FUNCTION:                gunelvemesh2VEMLab
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Read a gunelve mesh on a rectangular domain and write it into a VEMLab mesh format.
%
% Usage
% =====
% gunelvemesh2VEMLab(Node,Element,NElem,BoundaryNodes)
%
% Input
% =====
% Node    : array containing the nodal coordinates
% Element : cell array containing the element connectivity 
% NElem   : number of polygonal elements
% BoundaryNodes : structure containing the boundary nodes
% MeshFile : location and name of the file for writing the mesh
%
% Output
% ======
%
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Apr. 28, 2021: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gunelvemesh2VEMLab(Node,Element,NElem,BoundaryNodes,MeshFile)
  % for now, only rectangular domain is available
  gunelvemesh2VEMLab_rectangular_domain(Node,Element,NElem,BoundaryNodes,...
                                       MeshFile);
end

function gunelvemesh2VEMLab_rectangular_domain(Node,Element,NElem,...
                                              BoundaryNodes,MeshFile)
  fprintf('Printing mesh to a VEMLab mesh format...\n'); 
  fid = fopen(MeshFile,'w');
  % print domain type
  fprintf(fid,'# domain type\n');  
  fprintf(fid,'%s\n','RectangularDomain');    
  % print nodal coordinates
  fprintf(fid,'# nodal coordinates: number of nodes followed by the coordinates\n');
  nnode = size(Node,1);
  fprintf(fid,'%d\n',nnode);                                    
  for node_i = 1:nnode
    fprintf(fid,'%.16f %.16f\n', Node(node_i,1), Node(node_i,2));  
  end
  % print element connectivity
  fprintf(fid,'# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 node 3 ... node NVertex)\n');  
  fprintf(fid,'%d\n',NElem);                                 
  for el = 1:NElem
    NVertex = length(Element{el});
    fprintf(fid,'%d ', NVertex);
    for vertex = 1:(NVertex-1)
      fprintf(fid,'%d ', Element{el}(vertex));
    end
    fprintf(fid,'%d\n', Element{el}(NVertex));
  end
  % print bottom boundary  
  fprintf(fid,'# indices of nodes located on the bottom boundary\n');
  fprintf(fid,'%d ',BoundaryNodes.bottom);  
  fprintf(fid,'\n');
  % print top boundary  
  fprintf(fid,'# indices of nodes located on the top boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.top);
  fprintf(fid,'\n');  
  % print left boundary
  fprintf(fid,'# indices of nodes located on the left boundary\n');
  fprintf(fid,'%d ',BoundaryNodes.left);  
  fprintf(fid,'\n');  
  % print right boundary  
  fprintf(fid,'# indices of nodes located on the right boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.right);    
  fprintf(fid,'\n');
  % print xmin, xmax, ymin, ymax for the rectangular domain
  fprintf(fid,'# xmin, xmax, ymin, ymax of the bounding box\n'); 
  xmin=Node(BoundaryNodes.blcorner,1);
  xmax=Node(BoundaryNodes.trcorner,1);
  ymin=Node(BoundaryNodes.blcorner,2);
  ymax=Node(BoundaryNodes.trcorner,2); 
  fprintf(fid,'%.16f %.16f %.16f %.16f\n', xmin, xmax, ymin, ymax);   
  fclose(fid);
end

function plot_gunelve_mesh(connect,coords,config)

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
    patch('Faces',elements,'Vertices',coords,'FaceColor','w','LineWidth',config.plot_mesh_linewidth); 
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
    
    % plot nodes
    if strcmp(config.plot_mesh_nodes,'yes')
      hold on
      x=coords(:,1);
      y=coords(:,2);
      plot(x,y,'bo','MarkerFaceColor','b','MarkerSize',config.plot_mesh_nodesize);
      hold off
    end    

    % node numbering
    numnod=size(coords,1);
    if strcmp(config.plot_node_numbers,'yes')
      offset_x = 0.06*config.offset_hx; offset_y = 0.08*config.offset_hy;
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
