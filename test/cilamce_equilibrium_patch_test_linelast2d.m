%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         VEMLAB: A Matlab Library for the Virtual Element Method
%-------------------------------------------------------------------------------
%  Version      : 2.4                    
%  Source code  : http://camlab.cl/software/vemlab
%  Author       : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
%
%              (See Copyright and License notice in "license.txt")
%              (See updates and version details in "version.txt")
%-------------------------------------------------------------------------------
%
%               2D Linear Equilibrium Patch Test as Given in [2]
%                  Plane Strain - Linear Elastic Material
%
%-------------------------------------------------------------------------------
% References 
% ==========
% VEMLAB implementation is based on the VEM theory given in:
%
% [1] A. Ortiz-Bernardin, C. Alvarez, N. Hitschfeld-Kahler, A. Russo, 
%     R. Silva, A. Olate-Sanzana, "Veamy: an extensible object-oriented 
%     C++ library for the virtual element method," Numerical Algorithms, 
%     volume 82, pages 1189–1220(2019)
%
% [2] https://publicacoes.softaliza.com.br/cilamce2023/article/view/5017/4101
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  cilamce_equilibrium_patch_test_linelast2d

  close all;
  clear all;
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
  Ey=1e6;                               % Young's modulus
  nu=0.3;                                % Poisson's ratio
  elastic_state='plane_stress';          % 'plane_strain' or 'plane_stress'  
  
  % mesh filename: see available sample mesh files in folder "test/mesh_files/"
  %
  %     VEM2D meshes have the keyword "poly" in the mesh filename
  %     FEM2DQ4 meshes have the keyword "q4" in the mesh filename
  %     FEM2DT3 meshes have the keyword "t3" in the mesh filename 
  % 
  %     VEM2D meshes can be used with 'VEM2D' method only.
  %     FEM2DQ4 meshes can be used with 'FEM2DQ4' and 'VEM2D' methods.
  %     FEM2DT3 meshes can be used with 'FEM2DT3' and 'VEM2D' methods.
  %
  
  mesh_filename='cilamce_patch_test_q4_uniform_5x5elems.txt';  
  
  % method
  vemlab_method='VEM2D';        % 'VEM2D' (polygons - linear VEM) or
                                % 'FEM2DT3' (3-node triangles - linear FEM) or
                                % 'FEM2DQ4' (4-node quadrilaterals - bilinear FEM)
                                
  % Stability
  stability_type=2;           % 0 -> without stability
                              % 1 -> Gain et al.
                              % 2 -> D-recipe
                              % 3 -> modified D-recipe
                           
  % module
  vemlab_module='LinearElastostatics';  % 'LinearElastostatics' or 'Poisson'
  
  % solver
  vemlab_solver='sparse';   % 'sparse' or 'dense'    
  
  %% PLOT AND OUTPUT OPTIONS
  % to setup plot and output options go to folder 'config' and 
  % set the parameters inside function 'plot_and_output_options.m'
  
  %% VEMLAB CONFIGURATION
  
  config=config_vemlab(opsystem,vemlab_root_dir,mesh_filename,vemlab_module,...
                       vemlab_method,vemlab_solver,stability_type);
                                         
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
  
  body_force_fun_values.bx=@(x,y,matProps)bx_body_force_fun(x,y,matProps);
  body_force_fun_values.by=@(x,y,matProps)by_body_force_fun(x,y,matProps);
  
  %% ASSEMBLY 
  
  tic
  fprintf('Assemblying element matrices...\n'); 
  [K_global,f_global]=assembly(domainMesh,config,matProps,body_force_fun_values);
  
  %% NEUMANN BOUNDARY NODES/DOFS/FUNCTIONS (right side of the domain) 
  % (see definition of functions at the bottom of this file)

  Neumann_boundary_nodes=domainMesh.boundary_nodes.right;
  Neumann_boundary_dofs=domainMesh.boundary_dofs.right; 
  Neumann_fun_values.fx=@(x,y,matProps)fx_Neumann_fun(x,y,matProps);  
  Neumann_fun_values.fy=@(x,y,matProps)fy_Neumann_fun(x,y,matProps);   
  
  %% APPLY NEUMANN BCS ON THE BOUNDARY NODES 
  
  fprintf('Applying Neumann boundary conditions...\n');   
  Neumann_BCs=compute_Neumann_BCs(domainMesh,config,Neumann_boundary_nodes,...
                                  Neumann_boundary_dofs,Neumann_fun_values,matProps); 
  f_global(Neumann_BCs.indexes)=f_global(Neumann_BCs.indexes)+Neumann_BCs.values;
  
  %% DIRICHLET BOUNDARY NODES/DOFS/FUNCTIONS 
  % (see definition of functions at the bottom of this file)
  
  % left boundary
  Dirichlet_boundary_nodes1=domainMesh.boundary_nodes.left;
  Dirichlet_boundary_dofs1=domainMesh.boundary_dofs.left; 
  Dirichlet_fun_values1.ux=@(x,y,matProps)ux_Dirichlet_fun1(x,y,matProps);  
  Dirichlet_fun_values1.uy=@(x,y,matProps)uy_Dirichlet_fun1(x,y,matProps);   
  
  % bottom boundary
  Dirichlet_boundary_nodes2=domainMesh.boundary_nodes.bottom;
  Dirichlet_boundary_dofs2=domainMesh.boundary_dofs.bottom; 
  Dirichlet_fun_values2.ux=@(x,y,matProps)ux_Dirichlet_fun2(x,y,matProps);  
  Dirichlet_fun_values2.uy=@(x,y,matProps)uy_Dirichlet_fun2(x,y,matProps); 

  %% ENFORCE DIRICHLET BCS ON THE BOUNDARY NODES
  
  fprintf('Enforcing Dirichlet boundary conditions...\n');   
  u_nodal_sol=zeros(2*length(domainMesh.coords),1);   
  
  DB_dofs1=compute_Dirichlet_BCs(domainMesh,config,Dirichlet_boundary_nodes1,...
                                Dirichlet_boundary_dofs1,Dirichlet_fun_values1,matProps); % DOFs with Dirichlet BCs
  DB_dofs2=compute_Dirichlet_BCs(domainMesh,config,Dirichlet_boundary_nodes2,...
                                Dirichlet_boundary_dofs2,Dirichlet_fun_values2,matProps); % DOFs with Dirichlet BCs      
                              
  DB_dofs_indexes=[DB_dofs1.indexes;DB_dofs2.indexes];
  DB_dofs_values=[DB_dofs1.values;DB_dofs2.values];                                                      
  u_nodal_sol(DB_dofs_indexes)=DB_dofs_values;
  f_global=f_global-K_global*u_nodal_sol;    
  
  %% SOLVE RESULTING SYSTEM OF LINEAR EQUATIONS 
  
  fprintf('Solving system of linear equations...\n');   
  num_nodes=length(domainMesh.coords(:,1));    
  free_dofs=setdiff(1:2*num_nodes,DB_dofs_indexes); % dofs that do not have Dirichlet boundary conditions associated  
  u_nodal_sol(free_dofs)=K_global(free_dofs,free_dofs)\f_global(free_dofs); % nodal solution at free dofs
  toc
  
  %% EXACT SOLUTIONS FOR THE QUADRATIC PATCH TEST
  % (see definition of functions at the bottom of this file)
  
  exact_solution_handle.ux=@(x,y,matProps)ux_exact(x,y,matProps);
  exact_solution_handle.uy=@(x,y,matProps)uy_exact(x,y,matProps);
  exact_solution_handle.p=@(x,y,matProps)p_exact(x,y,matProps);   
  exact_solution_handle.strainvec=@(x,y,matProps)strainvec_exact(x,y,matProps);  
  
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

%% DEFINITION OF THE BODY FORCE FUNCTIONS

function bx = bx_body_force_fun(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the return "bx"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "bx" is an array that has the same form of the input "x" or "y"
  
%   N=length(x); % size of the b vector!
%   Ey=10^7;
%   nu=0.3;
%   D=(Ey/((1+nu)*(1-2*nu)))*[1-nu,nu,0; nu,1-nu,0; 0,0,2*(1-2*nu)];
%   bx=-0.2*D(1,1)-0.29*D(1,2)-0.4*D(3,3)*0.25; 
  bx=0;
end
function by = by_body_force_fun(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the return "by"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "by" is an array that has the same form of the input "x" or "y"
  % 
  by=0;  
end

%% DEFINITION OF NEUMANN FUNCTIONS

function fx = fx_Neumann_fun(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the return "fx"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "fx" is an array that has the same form of the input "x" or "y"
  
  fx=1000;  % is used as a force per length   
end
function fy = fy_Neumann_fun(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the return "fy"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "fy" is an array that has the same form of the input "x" or "y"
  
  fy=0;  % is used as a force per length   
end

%% DEFINITION OF DIRICHLET FUNCTIONS

function ux = ux_Dirichlet_fun1(x,y,matProps)
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
end
function uy = uy_Dirichlet_fun1(x,y,matProps)
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
  uy(:,2)=1; % free the y-dofs                       
end
function ux = ux_Dirichlet_fun2(x,y,matProps)
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
  ux(:,2)=1; % free the x-dofs                            
end
function uy = uy_Dirichlet_fun2(x,y,matProps)
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
end

%% DEFINITION OF THE EXACT SOLUTIONS

function ux = ux_exact(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "ux" is an array that has the same form of the input "x" or "y"
  Ey=matProps.Ey;
  nu=matProps.nu;
  fx = 1000;
  ux=(fx/Ey)*x;    
end
function uy = uy_exact(x,y,matProps)  
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "uy" is an array that has the same form of the input "x" or "y"
  
  Ey=matProps.Ey;
  nu=matProps.nu;
  fx = 1000;
  uy=-nu*fx/Ey*y;  
end
function p = p_exact(x,y,matProps)  
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "p" is an array that has the same form of the input "x" or "y"
    
  Ey=matProps.Ey;
  nu=matProps.nu;  
  D=matProps.D; % to be used with the strain in the format [e11 e22 e12]
  strain_vec=strainvec_exact(x,y,matProps);  
  stress_vec=D*strain_vec';
  if strcmp(matProps.plane_state,'plane_stress')
    p = -1/3*(stress_vec(1) + stress_vec(2));
  elseif strcmp(matProps.plane_state,'plane_strain')
    p = -1/3*(stress_vec(1) + stress_vec(2) + nu*(stress_vec(1) + stress_vec(2)));
  end

end
function strainvec = strainvec_exact(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "strainvec" is an array that has the same form of the input "x" or "y"
  
  Ey=matProps.Ey;
  nu=matProps.nu;
  fx=1000;
  duxdx=fx/Ey;  
  duydx=0;
  duxdy=0;
  duydy=-nu*fx/Ey;
 
  % strain vector in the form: [e11; e22; e12]  
  % note that the tird component is e12 and not 2*e12
  % this is due to that this code uses the material matrix D accordingly
  % with the strain vector definition as [e11; e22; e12]
  % -----------------
  % IMPORTANT: FOR x and y vectors of size n, the strain vector must be stored as: 
  % [e11_1 e22_1 e12_1;...;e11_n e22_n e12_n]
  % ------------------
  strainvec = [duxdx,duydy,0.5*(duxdy+duydx)]; 
end

