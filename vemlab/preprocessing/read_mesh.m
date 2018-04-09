%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                         read_mesh 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Read a mesh from a text file
%
% Usage
% =====
% domainMesh = read_mesh(config)
%
% Input
% =====
% config : structure containing program options and behavior
%
% Output
% ======
% domainMesh : a structure containing
%        domainMesh.coords         : nodal coordinates
%        domainMesh.connect        : connectivity of the elements in the mesh
%        domainMesh.boundary_nodes : nodes located on the boundary of the domain
%        domainMesh.boundary_dofs  : degrees of freedom of the boundary nodes
%
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
% Mar. 17, 2018: update to accomodate Poisson module (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

function domainMesh = read_mesh(config)
  fprintf('Reading a mesh...\n'); 
  meshFile = [config.mesh_folder_location,config.mesh_filename];
  % for now, only rectangular domain is available
  if strcmp(config.vemlab_module,'LinearElastostatics')
    domainMesh = read_mesh_rectangular_domain_linelast2d(meshFile);
  elseif strcmp(config.vemlab_module,'Poisson')
    domainMesh = read_mesh_rectangular_domain_poisson2d(meshFile);
  else
    throw_error('In read_mesh.m: vemlab_module\n');
  end
end

function domainMesh = read_mesh_rectangular_domain_linelast2d(meshFile)
  mesh_file = fopen(meshFile);
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of nodes that are located on the bottom boundary
  line = fgets(mesh_file); % read commented line
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.bottom = sscanf(line,'%d');
  % read indices of nodes that are located on the top boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.top = sscanf(line,'%d');
  % read indices of nodes that are located on the left boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.left = sscanf(line,'%d');    
  % read indices of nodes that are located on the right boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.right = sscanf(line,'%d');  
  % all boundary nodes
  domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.bottom;domainMesh.boundary_nodes.top;...
                             domainMesh.boundary_nodes.left;domainMesh.boundary_nodes.right];
  domainMesh.boundary_nodes.all=unique(domainMesh.boundary_nodes.all); % delete repeated nodes from boundary_nodes.all
  % read xmax, xmin, ymax, ymin for the rectangular domain
  line = fgets(mesh_file);
  domainMesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
  
  % form the boundary dofs
  % all
  domainMesh.boundary_dofs.all=zeros(2*length(domainMesh.boundary_nodes.all),1);      
  range=1:length(domainMesh.boundary_nodes.all);
  domainMesh.boundary_dofs.all(2*range-1)=2*domainMesh.boundary_nodes.all-1;
  domainMesh.boundary_dofs.all(2*range)=2*domainMesh.boundary_nodes.all;    
  % bottom
  domainMesh.boundary_dofs.bottom=zeros(2*length(domainMesh.boundary_nodes.bottom),1);      
  range=1:length(domainMesh.boundary_nodes.bottom);
  domainMesh.boundary_dofs.bottom(2*range-1)=2*domainMesh.boundary_nodes.bottom-1;
  domainMesh.boundary_dofs.bottom(2*range)=2*domainMesh.boundary_nodes.bottom;    
  % top
  domainMesh.boundary_dofs.top=zeros(2*length(domainMesh.boundary_nodes.top),1);      
  range=1:length(domainMesh.boundary_nodes.top);
  domainMesh.boundary_dofs.top(2*range-1)=2*domainMesh.boundary_nodes.top-1;
  domainMesh.boundary_dofs.top(2*range)=2*domainMesh.boundary_nodes.top;   
  % left
  domainMesh.boundary_dofs.left=zeros(2*length(domainMesh.boundary_nodes.left),1);      
  range=1:length(domainMesh.boundary_nodes.left);
  domainMesh.boundary_dofs.left(2*range-1)=2*domainMesh.boundary_nodes.left-1;
  domainMesh.boundary_dofs.left(2*range)=2*domainMesh.boundary_nodes.left;   
  % right
  domainMesh.boundary_dofs.right=zeros(2*length(domainMesh.boundary_nodes.right),1);      
  range=1:length(domainMesh.boundary_nodes.right);
  domainMesh.boundary_dofs.right(2*range-1)=2*domainMesh.boundary_nodes.right-1;
  domainMesh.boundary_dofs.right(2*range)=2*domainMesh.boundary_nodes.right;     
end

function domainMesh = read_mesh_rectangular_domain_poisson2d(meshFile)
  mesh_file = fopen(meshFile);
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    domainMesh.coords(i,1) = coordinates(1); domainMesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  domainMesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    domainMesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of nodes that are located on the bottom boundary
  line = fgets(mesh_file); % read commented line
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.bottom = sscanf(line,'%d');
  % read indices of nodes that are located on the top boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.top = sscanf(line,'%d');
  % read indices of nodes that are located on the left boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.left = sscanf(line,'%d');    
  % read indices of nodes that are located on the right boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  domainMesh.boundary_nodes.right = sscanf(line,'%d');  
  % all boundary nodes
  domainMesh.boundary_nodes.all = [domainMesh.boundary_nodes.bottom;domainMesh.boundary_nodes.top;...
                             domainMesh.boundary_nodes.left;domainMesh.boundary_nodes.right];
  domainMesh.boundary_nodes.all=unique(domainMesh.boundary_nodes.all); % delete repeated nodes from boundary_nodes.all
  % read xmax, xmin, ymax, ymin for the rectangular domain
  line = fgets(mesh_file);
  domainMesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
  
  % form the boundary dofs
  % all
  domainMesh.boundary_dofs.all=zeros(length(domainMesh.boundary_nodes.all),1);      
  range=1:length(domainMesh.boundary_nodes.all);
  domainMesh.boundary_dofs.all(range)=domainMesh.boundary_nodes.all;  
  % bottom
  domainMesh.boundary_dofs.bottom=zeros(length(domainMesh.boundary_nodes.bottom),1);      
  range=1:length(domainMesh.boundary_nodes.bottom);
  domainMesh.boundary_dofs.bottom(range)=domainMesh.boundary_nodes.bottom;   
  % top
  domainMesh.boundary_dofs.top=zeros(length(domainMesh.boundary_nodes.top),1);      
  range=1:length(domainMesh.boundary_nodes.top);
  domainMesh.boundary_dofs.top(range)=domainMesh.boundary_nodes.top;   
  % left
  domainMesh.boundary_dofs.left=zeros(length(domainMesh.boundary_nodes.left),1);      
  range=1:length(domainMesh.boundary_nodes.left);
  domainMesh.boundary_dofs.left(range)=domainMesh.boundary_nodes.left;  
  % right
  domainMesh.boundary_dofs.right=zeros(length(domainMesh.boundary_nodes.right),1);      
  range=1:length(domainMesh.boundary_nodes.right);
  domainMesh.boundary_dofs.right(range)=domainMesh.boundary_nodes.right;    
end