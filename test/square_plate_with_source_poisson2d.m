%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   VEMLab
%-------------------------------------------------------------------------------                                  
%  Version      : 2.0.2                        
%  Date         : May 13, 2018
%  Source code  : http://camlab.cl/research/software/vemlab
%  Author       : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
%
%              (See Copyright and License notice in "License.txt")
%              (See updates and version details in "Version.txt")
%-------------------------------------------------------------------------------
%
%                      Square Plate with Heat Source
%                           2D Poisson Problem
%
%-------------------------------------------------------------------------------
% References 
% ==========
% [1] A. Ortiz-Bernardin, C. Alvarez, N. Hitschfeld-Kahler, A. Russo, 
%     R. Silva, A. Olate-Sanzana, "Veamy: an extensible object-oriented 
%     C++ library for the virtual element method," arXiv:1708.03438 [cs.MS]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function square_plate_with_source_poisson2d

  close all;
  clc;
  
  %% THIS BLOCK MUST BE PUT ON EVERY TEST FILE
  
  % add all vemlab folders to the path
  opsystem=computer;
  is_Windows = strcmp(opsystem,'PCWIN') || strcmp(opsystem,'PCWIN64');
  is_Linux = strcmp(opsystem,'GLNX86') || strcmp(opsystem,'GLNXA64');   
  if is_Windows
    cd ..\; vemlab_root_dir=setpath;
  elseif is_Linux
    cd ../; vemlab_root_dir=setpath;
  end   
  
  %% INPUT DATA 
  
  % set material parameters
  k=1;                               % conductivity (isotropic material)
  
  % mesh filename: see available sample mesh files in folder "test/mesh_files/"
  %
  %     VEM2D meshes have the keyword "polygon" in the mesh filename
  %     FEM2DQ4 meshes have the keyword "q4" in the mesh filename
  %     FEM2DT3 meshes have the keyword "q4" in the mesh filename 
  % 
  %     VEM2D meshes can be used with 'VEM2D' method only.
  %     FEM2DQ4 meshes can be used with 'FEM2DQ4' and 'VEM2D' methods.
  %     FEM2DT3 meshes can be used with 'FEM2DT3' and 'VEM2D' methods.
  %
  mesh_filename='square_plate_poisson2d_1000poly_elems.txt';  

  % method
  vemlab_method='VEM2D';      % 'VEM2D' (polygons - linear VEM) or
                              % 'FEM2DT3' (3-node triangles - linear FEM) or
                              % 'FEM2DQ4' (4-node quadrilaterals - bilinear FEM)
  
  % module
  vemlab_module='Poisson';    % 'LinearElastostatics' or 'Poisson'
  
  %% PLOT AND OUTPUT OPTIONS
  % to setup plot and output options go to folder 'config' and 
  % set the parameters inside function 'plot_and_output_options.m'    
  
  %% VEMLAB CONFIGURATION
  
  config=config_vemlab(opsystem,vemlab_root_dir,mesh_filename,vemlab_module,...
                       vemlab_method);
                                         
  %% PRINT INIT MESSAGE TO SCREEN

  init_message(config);

  %% PREPROCESSING
  
  % read mesh
  domainMesh=read_mesh(config);  
  
  % plot mesh  
  plot_mesh2d(domainMesh,config);

  fprintf('\n');
  fprintf('Simulation started...\n');  
  
  %% CREATE POISSON MATERIAL 
  
  matProps=material_parameters_poisson(k);
  
  %% SOURCE TERM FUNCTIONS
  
  source_term_fun_values=@(x,y)source_term_fun(x,y);
  
  %% ASSEMBLY 
  
  tic
  fprintf('Assemblying element matrices...\n'); 
  [K_global,f_global]=assembly(domainMesh,config,matProps,source_term_fun_values);
  
  %% DIRICHLET BOUNDARY NODES/DOFS/FUNCTIONS
  % (see definition of functions at the bottom of this file)
  % The current problem has the entire boundary with Dirichlet BCs. Therefore,
  % we set the "Dirichlet_boundary_nodes" variable to "mesh.boundary_nodes.all".
  
  Dirichet_boundary_nodes=domainMesh.boundary_nodes.all; 
  Dirichet_boundary_dofs=domainMesh.boundary_dofs.all; 
  Dirichlet_fun_values=@(x,y)u_Dirichlet_fun(x,y);   

  %% ENFORCE DIRICHLET BCS ON THE BOUNDARY NODES
  
  fprintf('Enforcing Dirichlet boundary conditions...\n');   
  u_nodal_sol=zeros(length(domainMesh.coords),1);    
  DB_dofs=compute_Dirichlet_BCs(domainMesh,config,Dirichet_boundary_nodes,...
                                Dirichet_boundary_dofs,Dirichlet_fun_values); % DOFs with Dirichlet BCs
  u_nodal_sol(DB_dofs.indexes)=DB_dofs.values;
  f_global=f_global-K_global*u_nodal_sol;    
  
  %% SOLVE RESULTING SYSTEM OF LINEAR EQUATIONS 
  
  fprintf('Solving system of linear equations...\n');   
  num_nodes=length(domainMesh.coords(:,1));    
  free_dofs=setdiff(1:num_nodes,DB_dofs.indexes); % dofs that do not have Dirichlet boundary conditions associated  
  u_nodal_sol(free_dofs)=K_global(free_dofs,free_dofs)\f_global(free_dofs); % nodal solution at free dofs
  toc
  
  %% EXACT SOLUTIONS
  % (see definition of functions at the bottom of this file)
  
  exact_solution_handle.u=@(x,y)u_exact(x,y);
  exact_solution_handle.dudx=@(x,y)dudx_exact(x,y);   
  exact_solution_handle.dudy=@(x,y)dudy_exact(x,y);    
  
  %% NORMS OF THE ERROR
  
  fprintf('\n');
  fprintf('Computing norms of the solution error...\n');   
  compute_norms_of_the_error(exact_solution_handle,u_nodal_sol,domainMesh,config,matProps);  
  
  %% POSTPROCESSING
  
  % plot numerical solutions  
  postprocess_numerical_solution_poisson2d(domainMesh,u_nodal_sol,matProps,config);
  
  % plot exact solutions    
  postprocess_exact_solution_poisson2d(domainMesh,exact_solution_handle,matProps,config);  
  
  %% PRINT END MESSAGE TO SCREEN
  end_message;  
  
end

%% DEFINITION OF THE SOURCE TERM FUNCTION

function b = source_term_fun(x,y)
  % This function is intended to return a vector of size length(x) or length(y).
  % So, in case length(x) = 1, the return is just a single value.
  % Use something like x.*y (i.e., use the dot symbol), so that in case x and y 
  % are vectors, the return is a vector.
  b=32*y.*(1-y) + 32*x.*(1-x);
end


%% DEFINITION OF DIRICHLET FUNCTIONS FOR THE QUADRATIC PATCH TEST

function u = u_Dirichlet_fun(x,y)
  % x,y are vectors containing the coordinates of the nodes lying on the 
  % Dirichlet boundary, therefore this function is intended to return a vector
  % of size length(x) or length(y).
  % Use something like x.*y (i.e., use the dot symbol), so that in case x and y 
  % are vectors, the return is a vector.
  N=length(x);
  u=zeros(N,1);  
end

%% DEFINITION OF THE EXACT SOLUTIONS

function u = u_exact(x,y)
  % This function is intended to return a vector of size length(x) or length(y).
  % So, in case length(x) = 1, the return is just a single value.
  % Use something like x.*y (i.e., use the dot symbol), so that in case x and y 
  % are vectors, the return is a vector.  
  u=16*x.*y.*(1-x).*(1-y);
end

function dudx = dudx_exact(x,y)
  % This function is intended to return a vector of size length(x) or length(y).
  % So, in case length(x) = 1, the return is just a single value.
  % Use something like x.*y (i.e., use the dot symbol), so that in case x and y 
  % are vectors, the return is a vector. 
  dudx=16*y.*(1-y).*(1-2*x);
end

function dudy = dudy_exact(x,y)
  % This function is intended to return a vector of size length(x) or length(y).
  % So, in case length(x) = 1, the return is just a single value.
  % Use something like x.*y (i.e., use the dot symbol), so that in case x and y 
  % are vectors, the return is a vector.  
  dudy=16*x.*(1-x).*(1-2*y);
end



