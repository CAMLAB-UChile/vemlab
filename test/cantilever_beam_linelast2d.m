%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   VEMLab
%-------------------------------------------------------------------------------                                  
%  Version      : 2.2.2                      
%  Date         : 25-OCT-2019
%  Source code  : http://camlab.cl/software/vemlab
%  Author       : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
%
%              (See Copyright and License notice in "license.txt")
%              (See updates and version details in "version.txt")
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
  
  mesh_filename='cantilever_beam_250poly_elems.txt';   
  
  % method
  vemlab_method='VEM2D';      % 'VEM2D' (polygons - linear VEM) or
                              % 'FEM2DT3' (3-node triangles - linear FEM) or
                              % 'FEM2DQ4' (4-node quadrilaterals - bilinear FEM)
  
  % module
  vemlab_module='LinearElastostatics';  % 'LinearElastostatics' or 'Poisson'
  
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
  
  %% NEUMANN BOUNDARY NODES/DOFS/FUNCTIONS (right side of the beam) 
  % (see definition of functions at the bottom of this file)

  Neumann_boundary_nodes=domainMesh.boundary_nodes.right;
  Neumann_boundary_dofs=domainMesh.boundary_dofs.right; 
  Neumann_fun_values.fx=@(x,y)fx_Neumann_fun(x,y);  
  Neumann_fun_values.fy=@(x,y)fy_Neumann_fun(x,y);   
  
  %% APPLY NEUMANN BCS ON THE BOUNDARY NODES 
  
  fprintf('Applying Neumann boundary conditions...\n');   
  Neumann_BCs=compute_Neumann_BCs(domainMesh,config,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_values); 
  f_global(Neumann_BCs.indexes)=f_global(Neumann_BCs.indexes)+Neumann_BCs.values;
  
  %% DIRICHLET BOUNDARY NODES/DOFS/FUNCTIONS (left side of the beam) 
  % (see definition of functions at the bottom of this file)
  
  Dirichet_boundary_nodes=domainMesh.boundary_nodes.left;
  Dirichet_boundary_dofs=domainMesh.boundary_dofs.left; 
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
  
  %% EXACT SOLUTIONS FOR THE CANTILEVER BEAM
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

%% DEFINITION OF THE BODY FORCE FUNCTIONS FOR THE CANTILEVER BEAM

function bx = bx_body_force_fun(x,y)
  % Use something like x.*y (i.e., use the dot symbol) if the return "bx"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "bx" is an array that has the same form of the input "x" or "y"
  
  bx=0;  % is used as a force per volume 
end
function by = by_body_force_fun(x,y)
  % Use something like x.*y (i.e., use the dot symbol) if the return "by"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "by" is an array that has the same form of the input "x" or "y"
  
  by=0;  % is used as a force per volume  
end

%% DEFINITION OF NEUMANN FUNCTIONS FOR THE CANTILEVER BEAM

function fx = fx_Neumann_fun(x,y)
  % Use something like x.*y (i.e., use the dot symbol) if the return "fx"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "fx" is an array that has the same form of the input "x" or "y"
  
  fx=0;   % is used as a force per length 
end
function fy = fy_Neumann_fun(x,y)
  % Use something like x.*y (i.e., use the dot symbol) if the return "fy"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "fy" is an array that has the same form of the input "x" or "y"
  
  D=4; Ix=1*D*D*D/12; P=-1000;
  fy=P*(D*D/4 - y.*y)/(2*Ix);    % is used as a force per length 
end

%% DEFINITION OF DIRICHLET FUNCTIONS FOR THE CANTILEVER BEAM

function ux = ux_Dirichlet_fun(x,y)
  % INPUT: x,y are vectors containing the coordinates of the nodes lying on the 
  % Dirichlet boundary, therefore if the Dirichlet conditions depend on x and y,
  % consider using something like x.*y (i.e., use the dot symbol).
  %
  % OUTPUT: ux = array in which its first column contains the value of the 
  % Dirichlet boundary condition for the degrees of freedom-x, and its second 
  % column = free (1) or fixed (0) to indicate whether the corresponding value 
  % of the degree of freedom-x in the first column should be ignored (free) or 
  % applied (fixed).
 
  ux=zeros(length(x),2); % first column = value; second column = free (1) or fixed (0)
                         % In this case, all ux dofs are fixed and have the
                         % values of ux(:,1)
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  ux(:,1)=(-P/(6*Ey_bar*Ix))*y.*((6*L-3*x).*x +(2+nu_bar)*y.*y-3*D*D/2*(1+nu_bar));    
end
function uy = uy_Dirichlet_fun(x,y)
  % INPUT: x,y are vectors containing the coordinates of the nodes lying on the 
  % Dirichlet boundary, therefore if the Dirichlet conditions depend on x and y,
  % consider using something like x.*y (i.e., use the dot symbol).
  %
  % OUTPUT: uy = array in which its first column contains the value of the 
  % Dirichlet boundary condition for the degrees of freedom-y, and its second 
  % column = free (1) or fixed (0) to indicate whether the corresponding value 
  % of the degree of freedom-y in the first column should be ignored (free) or 
  % applied (fixed).
  
  uy=zeros(length(y),2); % first column = value; second column = free (1) or fixed (0)
                         % In this case, all uy dofs are fixed and have the
                         % values of uy(:,1)
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  uy(:,1)=P/(6*Ey_bar*Ix)*(3*nu_bar*y.*y.*(L-x)+(3*L-x).*x.*x);     
end

%% DEFINITION OF THE EXACT SOLUTIONS FOR THE CANTILEVER BEAM

function ux = ux_exact(x,y)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "ux" is an array that has the same form of the input "x" or "y"
  
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  ux=(-P/(6*Ey_bar*Ix))*y.*((6*L-3*x).*x +(2+nu_bar)*y.*y-3*D*D/2*(1+nu_bar));    
end
function uy = uy_exact(x,y)  
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "uy" is an array that has the same form of the input "x" or "y"
  
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  uy=P/(6*Ey_bar*Ix)*(3*nu_bar*y.*y.*(L-x)+(3*L-x).*x.*x);  
end
function duxdx = duxdx_exact(x,y)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "duxdx" is an array that has the same form of the input "x" or "y"
  
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  duxdx=-(P*y.*(6*L-6*x))/(6*Ey_bar*Ix);
end
function duydx = duydx_exact(x,y)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "duydx" is an array that has the same form of the input "x" or "y"
  
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  duydx=-(P*(3*nu_bar*y.^2-2*x.*(3*L-x)+x.^2))/(6*Ey_bar*Ix);
end
function duxdy = duxdy_exact(x,y)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "duxdy" is an array that has the same form of the input "x" or "y"
  
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  duxdy=-(P*((nu_bar+2)*y.^2+x.*(6*L-3*x)-(3*D^2*(nu_bar+1))/2))/(6*Ey_bar*Ix)...
        -(P*y.^2*(nu_bar+2))/(3*Ey_bar*Ix);
end
function duydy = duydy_exact(x,y)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "duydy" is an array that has the same form of the input "x" or "y"
  
  Ey=10^7; nu=0.3; Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  duydy=(P*nu_bar*y.*(L-x))/(Ey_bar*Ix);
end


