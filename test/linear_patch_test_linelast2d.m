%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   VEMLab
%-------------------------------------------------------------------------------                                  
%  Version      : 2.0.1                         
%  Date         : April 26, 2018
%  Source code  : http://camlab.cl/research/software/vemlab
%  Author       : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
%
%              (See Copyright and License notice in "License.txt")
%              (See updates and version details in "Version.txt")
%-------------------------------------------------------------------------------
%
%                           2D Linear Patch Test
%                  Plane Strain - Linear Elastic Material
%
%-------------------------------------------------------------------------------
% References 
% ==========
% [1] A. Ortiz-Bernardin, C. Alvarez, N. Hitschfeld-Kahler, A. Russo, 
%     R. Silva, A. Olate-Sanzana, "Veamy: an extensible object-oriented 
%     C++ library for the virtual element method," arXiv:1708.03438 [cs.MS]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function linear_patch_test_linelast2d

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
  Ey=10^7;                               % Young's modulus
  nu=0.3;                                % Poisson's ratio
  elastic_state='plane_strain';          % 'plane_strain' or 'plane_stress'  
  
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
  mesh_filename='patch_test_50poly_elems.txt';  
  
  % method
  vemlab_method='VEM2D';      % 'VEM2D' (polygons - linear VEM) or
                              % 'FEM2DT3' (3-node triangles - linear FEM) or
                              % 'FEM2DQ4' (4-node quadrilaterals - bilinear FEM)
  
  % module
  vemlab_module='LinearElastostatics';   % 'LinearElastostatics' or 'Poisson'
  
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
  
  %% CREATE LINEAR ELASTIC MATERIAL 
  
  matProps=material_parameters_linelast(Ey,nu,elastic_state);
  
  %% BODY FORCE FUNCTIONS
  
  body_force_fun_values.bx=@(x,y)bx_body_force_fun(x,y);
  body_force_fun_values.by=@(x,y)by_body_force_fun(x,y);
  
  %% ASSEMBLY 
  
  tic
  fprintf('Assemblying element matrices...\n'); 
  [K_global,f_global]=assembly(domainMesh,config,matProps,body_force_fun_values);
  
  
  %% DIRICHLET BOUNDARY NODES/DOFS/FUNCTIONS (entire boundary) 
  % (see definition of functions at the bottom of this file)
  
  Dirichet_boundary_nodes=domainMesh.boundary_nodes.all;
  Dirichet_boundary_dofs=domainMesh.boundary_dofs.all; 
  Dirichlet_fun_values.ux=@(x,y)ux_Dirichlet_fun(x,y);  
  Dirichlet_fun_values.uy=@(x,y)uy_Dirichlet_fun(x,y);   

  %% ENFORCE DIRICHLET BCS ON THE BOUNDARY NODES
  
  fprintf('Enforcing Dirichlet boundary conditions...\n');   
  u_nodal_sol=zeros(2*length(domainMesh.coords),1);    
  DB_dofs=compute_Dirichlet_BCs(domainMesh,config,Dirichet_boundary_nodes,...
                                Dirichet_boundary_dofs,Dirichlet_fun_values); % DOFs with Dirichlet BCs
  u_nodal_sol(DB_dofs.indexes)=DB_dofs.values;
  f_global=f_global-K_global*u_nodal_sol;    
  
  %% SOLVE RESULTING SYSTEM OF LINEAR EQUATIONS 
  
  fprintf('Solving system of linear equations...\n');   
  num_nodes=length(domainMesh.coords(:,1));    
  free_dofs=setdiff(1:2*num_nodes,DB_dofs.indexes); % dofs that do not have Dirichlet boundary conditions associated  
  u_nodal_sol(free_dofs)=K_global(free_dofs,free_dofs)\f_global(free_dofs); % nodal solution at free dofs
  toc
  
  %% EXACT SOLUTIONS FOR THE LINEAR PATCH TEST
  % (see definition of functions at the bottom of this file)
  
  exact_solution_handle.ux=@(x,y)ux_exact(x,y);
  exact_solution_handle.uy=@(x,y)uy_exact(x,y);
  exact_solution_handle.duxdx=@(x,y)duxdx_exact(x,y);  
  exact_solution_handle.duydx=@(x,y)duydx_exact(x,y);  
  exact_solution_handle.duxdy=@(x,y)duxdy_exact(x,y);  
  exact_solution_handle.duydy=@(x,y)duydy_exact(x,y);   
  
  %% NORMS OF THE ERROR
  
  fprintf('\n');
  fprintf('Computing norms of the solution error...\n');   
  compute_norms_of_the_error(exact_solution_handle,u_nodal_sol,domainMesh,config,matProps);  
  
  %% POSTPROCESSING
  
  % plot numerical solutions  
  postprocess_numerical_solution_linelast2d(domainMesh,u_nodal_sol,matProps,config);
  
  % plot exact solutions    
  postprocess_exact_solution_linelast2d(domainMesh,exact_solution_handle,matProps,config); 
  
  %% PRINT END MESSAGE TO SCREEN
  end_message;  
  
end

%% DEFINITION OF THE BODY FORCE FUNCTIONS FOR THE LINEAR PATCH TEST

function bx = bx_body_force_fun(x,y)
  bx=0;    
end
function by = by_body_force_fun(x,y)
  by=0;    
end

%% DEFINITION OF DIRICHLET FUNCTIONS FOR THE LINEAR PATCH TEST

function ux = ux_Dirichlet_fun(x,y)
  ux=x;    
end
function uy = uy_Dirichlet_fun(x,y)
  uy=x+y;     
end

%% DEFINITION OF THE EXACT SOLUTIONS FOR THE LINEAR PATCH TEST

function ux = ux_exact(x,y)
  ux=x;    
end
function uy = uy_exact(x,y)  
  uy=x+y;  
end
function duxdx = duxdx_exact(x,y)
  duxdx=1;
end
function duydx = duydx_exact(x,y)
  duydx=1;
end
function duxdy = duxdy_exact(x,y)
  duxdy=0;
end
function duydy = duydy_exact(x,y)
  duydy=1;
end


