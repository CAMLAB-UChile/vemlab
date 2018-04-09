%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                   read_mesh 
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
% [coords,connect,boundary_nodes] = read_mesh(meshFile)
%
% Input
% =====
% meshFile : name of the file containing the mesh (string)
%
% Output
% ======
% coords         : nodal coordinates
% connect        : connectivity of the elements in the mesh
% boundary_nodes : nodes located on the boundary of the domain
%
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

function mesh = read_mesh(meshFile)
  fprintf('Reading a mesh...\n'); 
  % for now, only rectangular domain is available
  mesh = read_mesh_rectangular_domain(meshFile);
end

function mesh = read_mesh_rectangular_domain(meshFile)
  mesh_file = fopen(meshFile);
  % read nodal coordinates
  line = fgets(mesh_file); % read commented line
  nodes_number = sscanf(fgets(mesh_file),'%d');
  mesh.coords = zeros(nodes_number,2);
  for i=1:nodes_number
    line = fgets(mesh_file);
    coordinates = sscanf(line,'%f');
    mesh.coords(i,1) = coordinates(1); mesh.coords(i,2) = coordinates(2);
  end
  % read connectivities
  line = fgets(mesh_file); % read commented line
  polygons_number = sscanf(fgets(mesh_file),'%d');
  mesh.connect = cell(polygons_number,1);
  for i=1:polygons_number
    line = fgets(mesh_file); 
    indexes = sscanf(line,'%d');
    mesh.connect{i,1} = indexes(2:indexes(1)+1);
  end
  % read indices of nodes that are located on the bottom boundary
  line = fgets(mesh_file); % read commented line
  line = fgets(mesh_file);
  mesh.boundary_nodes.bottom = sscanf(line,'%d');
  % read indices of nodes that are located on the top boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  mesh.boundary_nodes.top = sscanf(line,'%d');
  % read indices of nodes that are located on the left boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  mesh.boundary_nodes.left = sscanf(line,'%d');    
  % read indices of nodes that are located on the right boundary  
  line = fgets(mesh_file); % read commented line  
  line = fgets(mesh_file);
  mesh.boundary_nodes.right = sscanf(line,'%d');  
  % all boundary nodes
  mesh.boundary_nodes.all = [mesh.boundary_nodes.bottom;mesh.boundary_nodes.top;...
                             mesh.boundary_nodes.left;mesh.boundary_nodes.right];
  mesh.boundary_nodes.all=unique(mesh.boundary_nodes.all); % delete repeated nodes from boundary_nodes.all
  % read xmax, xmin, ymax, ymin for the rectangular domain
  line = fgets(mesh_file);
  mesh.BdBox = sscanf(line,'%f');  
  fclose(mesh_file);
  
  % form the boundary dofs
  % all
  mesh.boundary_dofs.all=zeros(2*length(mesh.boundary_nodes.all),1);      
  range=1:length(mesh.boundary_nodes.all);
  mesh.boundary_dofs.all(2*range-1)=2*mesh.boundary_nodes.all-1;
  mesh.boundary_dofs.all(2*range)=2*mesh.boundary_nodes.all;    
  % bottom
  mesh.boundary_dofs.bottom=zeros(2*length(mesh.boundary_nodes.bottom),1);      
  range=1:length(mesh.boundary_nodes.bottom);
  mesh.boundary_dofs.bottom(2*range-1)=2*mesh.boundary_nodes.bottom-1;
  mesh.boundary_dofs.bottom(2*range)=2*mesh.boundary_nodes.bottom;    
  % top
  mesh.boundary_dofs.top=zeros(2*length(mesh.boundary_nodes.top),1);      
  range=1:length(mesh.boundary_nodes.top);
  mesh.boundary_dofs.top(2*range-1)=2*mesh.boundary_nodes.top-1;
  mesh.boundary_dofs.top(2*range)=2*mesh.boundary_nodes.top;   
  % left
  mesh.boundary_dofs.left=zeros(2*length(mesh.boundary_nodes.left),1);      
  range=1:length(mesh.boundary_nodes.left);
  mesh.boundary_dofs.left(2*range-1)=2*mesh.boundary_nodes.left-1;
  mesh.boundary_dofs.left(2*range)=2*mesh.boundary_nodes.left;   
  % right
  mesh.boundary_dofs.right=zeros(2*length(mesh.boundary_nodes.right),1);      
  range=1:length(mesh.boundary_nodes.right);
  mesh.boundary_dofs.right(2*range-1)=2*mesh.boundary_nodes.right-1;
  mesh.boundary_dofs.right(2*range)=2*mesh.boundary_nodes.right;     
  
end