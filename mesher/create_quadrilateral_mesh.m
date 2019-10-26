function create_quadrilateral_mesh
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
  
  xmin=0; xmax=1; 
  ymin=0; ymax=1;
  ndiv_x=5; ndiv_y=5; % number of divisions along x and y coordinates
  mesh_type='uniform';  % 'uniform'
  mesh_filename='patch_test_q4_uniform_5x5elems.txt'; 
  
  %%%%%%%%%%%%%%%%%%%        END USER INPUT DATA       %%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%
  % configure mesher
  config=config_vemlab_mesher(opsystem,vemlab_root_dir,mesh_filename);
  
  % run quad4_mesh
  [coords,connect]=quad4mesh(mesh_type,xmin,ymin,xmax,ymax,ndiv_x,ndiv_y);
  
  % get boundary nodes
  BdBox=[xmin,ymin;xmax,ymax];   
  BoundaryNodes=BoundaryRectangularDomain(coords,BdBox);
  
  % distmesh2VEMLab: write mesh to a VEMLab mesh format
  mesh_file=[config.mesh_folder_location,mesh_filename];  
  quad4mesh2VEMLab(coords,connect,size(connect,1),BoundaryNodes,mesh_file);   
 
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
% FUNCTION:                quad4mesh2VEMLab
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Read a quad4mesh mesh and write it into a VEMLab mesh format.
%
% Usage
% =====
% quad4mesh2VEMLab(Node,Element,NElem,BoundaryNodes)
%
% Input
% =====
% Node    : quad4mesh array containing the nodal coordinates
% Element : quad4mesh cell array containing the element connectivity 
% NElem   : number of 4-node quadrilateral elements
% BoundaryNodes : Structure containing the boundary nodes
% MeshFile : Location and name of the file for writing the mesh
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
% May 17, 2018: add domainType string (by A. Ortiz-Bernardin)
% Jan. 8, 2018: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function quad4mesh2VEMLab(Node,Element,NElem,BoundaryNodes,MeshFile)
  % for now, only rectangular domain is available
  quad4mesh2VEMLab_rectangular_domain(Node,Element,NElem,BoundaryNodes,...
                                       MeshFile);
end

function quad4mesh2VEMLab_rectangular_domain(Node,Element,NElem,...
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
  fprintf(fid,'# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 node 3 node 4)\n');  
  fprintf(fid,'%d\n',NElem);                                 
  for el = 1:NElem
    NVertex = length(Element(el,:));
    fprintf(fid,'%d ', NVertex);
    for vertex = 1:(NVertex-1)
      fprintf(fid,'%d ', Element(el,vertex));
    end
    fprintf(fid,'%d\n', Element(el,NVertex));
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
