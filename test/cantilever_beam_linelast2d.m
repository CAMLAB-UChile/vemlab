%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   VEMLab
%-------------------------------------------------------------------------------                                  
%  Version      : 1.0                         
%  Date         : 17-FEB-2018
%  Source code  : http://camlab.cl/research/software/vemlab
%  Author       : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
%
%              (See Copyright and License notice in "License.txt")
%              (See updates and version details in "Version.txt")
%-------------------------------------------------------------------------------
%
%             2D Cantilever Beam Subjected to a Parabolic End Load
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

function cantilever_beam_linelast2d

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
  
  %% INPUT DATA 
  
  % set material parameters
  Ey=10^7;                               % Young's modulus
  nu=0.3;                                % Poisson's ratio
  elastic_state='plane_strain';          % 'plane_strain' or 'plane_stress'  
  
  % mesh filename 
  mesh_filename='cantilever_beam_50poly_elems.txt';  
  
  % options
  plot_mesh_over_results='no';           % 'yes' or 'no'
  
  % method
  vemlab_method='VEM2D';                 % 'VEM2D', 'FEM2DT3' or 'FEM2DQ4'
  
  % module
  vemlab_module='LinearElastostatics';  % 'LinearElastostatics'
  
  %% VEMLAB CONFIGURATION
  
  config=config_vemlab(opsystem,vemlab_root_dir,mesh_filename,...
                       plot_mesh_over_results,vemlab_module,vemlab_method);
                                         
  %% PRINT INIT MESSAGE TO SCREEN

  init_message(config);

  %% PREPROCESSING
  
  % read mesh
  mesh=read_mesh([config.mesh_folder_location,config.mesh_filename]);  
  
  % plot mesh  
  plot_mesh2d(mesh);

  fprintf('\n');
  fprintf('Simulation started...\n');  
  
  %% CREATE LINEAR ELASTIC MATERIAL 
  
  matProps=material_parameters_linelast(Ey,nu,elastic_state);
  
  %% BODY FORCE FUNCTIONS
  
  body_force_fun_value.bx=@(x,y)bx_body_force_fun(x,y);
  body_force_fun_value.by=@(x,y)by_body_force_fun(x,y);
  
  %% ASSEMBLY 
  
  tic
  fprintf('Assemblying element matrices...\n'); 
  [K_global,f_global]=assembly(mesh,config,matProps,body_force_fun_value);
  
  %% NEUMANN BOUNDARY NODES/DOFS/FUNCTIONS (right side of the beam) 
  % (see definition of functions at the bottom of this file)

  Neumann_boundary_nodes=mesh.boundary_nodes.right;
  Neumann_boundary_dofs=mesh.boundary_dofs.right; 
  Neumann_fun_value.fx=@(x,y)fx_Neumann_fun(x,y);  
  Neumann_fun_value.fy=@(x,y)fy_Neumann_fun(x,y);   
  
  %% APPLY NEUMANN BCS ON THE BOUNDARY NODES 
  
  fprintf('Applying Neumann boundary conditions...\n');   
  Neumann_BCs=compute_Neumann_BCs(mesh,config,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_value); 
  f_global(Neumann_BCs.indexes)=f_global(Neumann_BCs.indexes)+Neumann_BCs.values;
  
  %% DIRICHLET BOUNDARY NODES/DOFS/FUNCTIONS (left side of the beam) 
  % (see definition of functions at the bottom of this file)
  
  Dirichet_boundary_nodes=mesh.boundary_nodes.left;
  Dirichet_boundary_dofs=mesh.boundary_dofs.left; 
  Dirichlet_fun_value.ux=@(x,y)ux_Dirichlet_fun(x,y);  
  Dirichlet_fun_value.uy=@(x,y)uy_Dirichlet_fun(x,y);   

  %% ENFORCE DIRICHLET BCS ON THE BOUNDARY NODES
  
  fprintf('Enforcing Dirichlet boundary conditions...\n');   
  u_nodal_sol=zeros(2*length(mesh.coords),1);    
  DB_dofs=compute_Dirichlet_BCs(mesh,config,Dirichet_boundary_nodes,...
                                Dirichet_boundary_dofs,Dirichlet_fun_value); % DOFs with Dirichlet BCs
  u_nodal_sol(DB_dofs.indexes)=DB_dofs.values;
  f_global=f_global-K_global*u_nodal_sol;    
  
  %% SOLVE RESULTING SYSTEM OF LINEAR EQUATIONS 
  
  fprintf('Solving system of linear equations...\n');   
  num_nodes=length(mesh.coords(:,1));    
  free_dofs=setdiff(1:2*num_nodes,DB_dofs.indexes); % dofs that do not have Dirichlet boundary conditions associated  
  u_nodal_sol(free_dofs)=K_global(free_dofs,free_dofs)\f_global(free_dofs); % nodal solution at free dofs
  toc
  
  %% EXACT SOLUTIONS FOR THE CANTILEVER BEAM
  % (see definition of functions at the bottom of this file)
  
  exact_sol.ux=@(x,y)ux_exact(x,y);
  exact_sol.uy=@(x,y)uy_exact(x,y);
  exact_sol.duxdx=@(x,y)duxdx_exact(x,y);  
  exact_sol.duydx=@(x,y)duydx_exact(x,y);  
  exact_sol.duxdy=@(x,y)duxdy_exact(x,y);  
  exact_sol.duydy=@(x,y)duydy_exact(x,y);   
  
  %% NORMS OF THE ERROR
  
  fprintf('\n');
  fprintf('Computing norms of the solution error...\n');   
  compute_norms_of_the_error(exact_sol,u_nodal_sol,mesh,config,matProps);  
  
  %% POSTPROCESSING
  
  % compute stresses and strains
  [stresses,strains]=compute_stresses_and_strains(u_nodal_sol,mesh,matProps,config);
  
  % plot VEM solution  
  plot_solution_linelast2d(mesh,u_nodal_sol,'$u_x^{\textrm{vem}}$',...
                        '$u_y^{\textrm{vem}}$','$||u^{\textrm{vem}}||$',...
                        config.plot_mesh_over_results,config.vemlab_method);  
  % plot exact solution  
  [u_nodal_exact,~,~]=exact_solutions_linelast2d(exact_sol,mesh.coords);
  plot_solution_linelast2d(mesh,u_nodal_exact,'$u_x^{\textrm{exact}}$',...
                        '$u_y^{\textrm{exact}}$','$||u^{\textrm{exact}}||$',...
                        config.plot_mesh_over_results,'exact');
  
  % write VEM solution to a text file
  write_solution_txt_linelast2d(mesh,u_nodal_sol,stresses,strains,config);
  
  % write VEM solution to a GiD file
  write_solution_GiD_linelast2d(mesh,u_nodal_sol,stresses,strains,config) 
  
  %% PRINT END MESSAGE TO SCREEN
  end_message;  
  
end

%% DEFINITION OF THE BODY FORCE FUNCTIONS FOR THE CANTILEVER BEAM

function bx = bx_body_force_fun(x,y)
  bx=0;    
end
function by = by_body_force_fun(x,y)
  by=0;    
end

%% DEFINITION OF NEUMANN FUNCTIONS FOR THE CANTILEVER BEAM

function fx = fx_Neumann_fun(x,y)
  fx=0;    
end
function fy = fy_Neumann_fun(x,y)
  D=4; Ix=1*D*D*D/12; P=-1000;
  fy=P*(D*D/4 - y.*y)/(2*Ix);    
end

%% DEFINITION OF DIRICHLET FUNCTIONS FOR THE CANTILEVER BEAM

function ux = ux_Dirichlet_fun(x,y)
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  ux=(-P/(6*Ey_bar*Ix))*y.*((6*L-3*x).*x +(2+nu_bar)*y.*y-3*D*D/2*(1+nu_bar));    
end
function uy = uy_Dirichlet_fun(x,y)
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  uy=P/(6*Ey_bar*Ix)*(3*nu_bar*y.*y.*(L-x)+(3*L-x).*x.*x);     
end

%% DEFINITION OF THE EXACT SOLUTIONS FOR THE CANTILEVER BEAM

function ux = ux_exact(x,y)
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  ux=(-P/(6*Ey_bar*Ix))*y.*((6*L-3*x).*x +(2+nu_bar)*y.*y-3*D*D/2*(1+nu_bar));    
end
function uy = uy_exact(x,y)  
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  uy=P/(6*Ey_bar*Ix)*(3*nu_bar*y.*y.*(L-x)+(3*L-x).*x.*x);  
end
function duxdx = duxdx_exact(x,y)
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  duxdx=-(P*y.*(6*L-6*x))/(6*Ey_bar*Ix);
end
function duydx = duydx_exact(x,y)
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  duydx=-(P*(3*nu_bar*y.^2-2*x.*(3*L-x)+x.^2))/(6*Ey_bar*Ix);
end
function duxdy = duxdy_exact(x,y)
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  duxdy=-(P*((nu_bar+2)*y.^2+x.*(6*L-3*x)-(3*D^2*(nu_bar+1))/2))/(6*Ey_bar*Ix)...
        -(P*y.^2*(nu_bar+2))/(3*Ey_bar*Ix);
end
function duydy = duydy_exact(x,y)
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  duydy=(P*nu_bar*y.*(L-x))/(Ey_bar*Ix);
end


