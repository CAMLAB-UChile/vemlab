%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            VEMLab (extended version)
%                           Source code  : not released
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                         read_T3_mesh 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Read a mesh of 3-node triangular elements from a text file to be used for
% conversion to a polygonal mesh.
%
% Usage
% =====
% domainMesh = read_T3_mesh(config)
%
% Input
% =====
% config : structure containing program options and behavior
%
% Output
% ======
% domainMesh : a structure containing
%        domainMesh.domain_type    : domain type
%        domainMesh.num_nodes      : number of nodes in the T3 mesh
%        domainMesh.coords         : nodal coordinates
%        domainMesh.num_elements   : number of T3 elements
%        domainMesh.connect        : connectivity of the elements in the mesh
%        domainMesh.num_corner_nodes : number of corner nodes
%        domainMesh.corner_nodes   : nodes located on the corners of the domain
%        domainMesh.num_boundary_edges : number of boundary edges
%        domainMesh.boundary_edges : connectivity of the domain boundary edges            
%        domainMesh.BdBox          : domain bounding box [xmin,xmax,ymin,ymax]
%        domainMesh.num_dirichlet_entities : number of entities that belong
%                                            to the Dirichlet boundary
%
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% May 4, 2021: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

function domainMesh = read_T3_mesh(meshFile)
  mesh_file = fopen(meshFile);
  % read domain type 
  line = fgets(mesh_file); % read commented line  
  domain_type = sscanf(fgets(mesh_file),'%s');
  domainMesh.domain_type=domain_type;
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.num_nodes=nodes_number; 
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.num_elements=polygons_number;  
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of corner nodes: first entry is the node and second entry is an integer 
  % (0,1,2) that indicates whether the corner node is also a Dirichlet or Neumann node, 
  % i.e., 0=None, 1=Dirichlet, 2=Neumann, and the third entry is an integer that can be 
  % used as a tag to distinguish when 2 or more corner nodes will have different 
  % Dirichlet or Neumann boundary conditions
  domainMesh.num_dirichlet_entities = 0;
  domainMesh.num_neumann_entities = 0;  
  line = fgets(mesh_file); % read commented line
  domainMesh.num_corner_nodes=sscanf(fgets(mesh_file),'%d');
  domainMesh.corner_nodes=zeros(domainMesh.num_corner_nodes,3);
  for i=1:domainMesh.num_corner_nodes
    line = fgets(mesh_file);      
    domainMesh.corner_nodes(i,:) = sscanf(line,'%d')';
    if domainMesh.corner_nodes(i,2)==1 % Dirichlet boundary node
      domainMesh.num_dirichlet_entities = domainMesh.num_dirichlet_entities + 1;
    elseif domainMesh.corner_nodes(i,2)==2 % Neumann boundary node
      domainMesh.num_neumann_entities = domainMesh.num_neumann_entities + 1;     
    end
  end
  % read connectivity of boundary edges: first and second entries are the nodal connectivity;
  % the third entry is an integer (0,1,2) that indicates whether the edge is also a 
  % Dirichlet or Neumann boundary edge, i.e., 0=None, 1=Dirichlet, 2=Neumann, and the fourth 
  % entry is an integer that can be used as a tag to distinguish when 2 or more boundary 
  % edges will have different Dirichlet or Neumann boundary conditions
  line = fgets(mesh_file); % read commented line
  domainMesh.num_boundary_edges=sscanf(fgets(mesh_file),'%d');
  domainMesh.boundary_edges=zeros(domainMesh.num_boundary_edges,4);
  for i=1:domainMesh.num_boundary_edges
    line = fgets(mesh_file);
    domainMesh.boundary_edges(i,:) = sscanf(line,'%d')';
    if domainMesh.boundary_edges(i,3)==1 % Dirichlet boundary edge
      domainMesh.num_dirichlet_entities = domainMesh.num_dirichlet_entities + 1;
    elseif domainMesh.boundary_edges(i,3)==2 % Neumann boundary edge
      domainMesh.num_neumann_entities = domainMesh.num_neumann_entities + 1;      
    end    
  end   
  % calculate bounding box
  domainMesh.BdBox=...
    [min(domainMesh.coords(:,1)), max(domainMesh.coords(:,1)),...
    min(domainMesh.coords(:,2)), max(domainMesh.coords(:,2))];
end
