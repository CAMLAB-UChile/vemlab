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
%             2D Cantilever Beam Subjected to a Parabolic End Load
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cantilever_beam_linelast2d

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
  Ey=10^7;                               % Young's modulus
  nu=0.3;                                % Poisson's ratio
  elastic_state='plane_strain';          % 'plane_strain' or 'plane_stress'  
  
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
  
  mesh_filename='cantilever_beam_1500poly_elems.txt';   
%   mesh_filename='cantilever_beam_q4_uniform_66x33.txt';  
%   mesh_filename='cantilever_beam_t3_hsize_01.txt';    
  
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
  
  %% NEUMANN BOUNDARY NODES/DOFS/FUNCTIONS (right side of the beam) 
  % (see definition of functions at the bottom of this file)

  Neumann_boundary_nodes=domainMesh.boundary_nodes.right;
  Neumann_boundary_dofs=domainMesh.boundary_dofs.right; 
  Neumann_fun_values.fx=@(x,y,matProps)fx_Neumann_fun(x,y,matProps);  
  Neumann_fun_values.fy=@(x,y,matProps)fy_Neumann_fun(x,y,matProps);   
  
  %% APPLY NEUMANN BCS ON THE BOUNDARY NODES 
  
  fprintf('Applying Neumann boundary conditions...\n');   
  Neumann_BCs=compute_Neumann_BCs(domainMesh,config,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_values,matProps); 
  f_global(Neumann_BCs.indexes)=f_global(Neumann_BCs.indexes)+Neumann_BCs.values;
  
  %% DIRICHLET BOUNDARY NODES/DOFS/FUNCTIONS (left side of the beam) 
  % (see definition of functions at the bottom of this file)
  
  Dirichet_boundary_nodes=domainMesh.boundary_nodes.left;
  Dirichet_boundary_dofs=domainMesh.boundary_dofs.left; 
  Dirichlet_fun_values.ux=@(x,y,matProps)ux_Dirichlet_fun(x,y,matProps);  
  Dirichlet_fun_values.uy=@(x,y,matProps)uy_Dirichlet_fun(x,y,matProps);   

  %% ENFORCE DIRICHLET BCS ON THE BOUNDARY NODES
  
  fprintf('Enforcing Dirichlet boundary conditions...\n');   
  u_nodal_sol=zeros(2*length(domainMesh.coords),1);    
  DB_dofs=compute_Dirichlet_BCs(domainMesh,config,Dirichet_boundary_nodes,...
                                Dirichet_boundary_dofs,Dirichlet_fun_values,matProps); % DOFs with Dirichlet BCs
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
  
  exact_solution_handle.ux=@(x,y,matProps)ux_exact(x,y,matProps);
  exact_solution_handle.uy=@(x,y,matProps)uy_exact(x,y,matProps);
  exact_solution_handle.p=@(x,y,matProps)p_exact(x,y,matProps);   
  exact_solution_handle.strainvec=@(x,y,matProps)strainvec_exact(x,y,matProps);   
  
  %% NORMS OF THE ERROR
  
  fprintf('\n');
  fprintf('Computing norms of the solution error...\n');   
  [h_max,L2rel,H1rel,L2prel]=...
    compute_norms_of_the_error(exact_solution_handle,u_nodal_sol,...
                               domainMesh,config,matProps);  
  
  %% POSTPROCESSING
  
  % plot numerical solutions  
  postprocess_numerical_solution_linelast2d(domainMesh,u_nodal_sol,matProps,config);
  
  % plot exact solutions    
  postprocess_exact_solution_linelast2d(domainMesh,exact_solution_handle,matProps,config); 
  
  %% PRINT END MESSAGE TO SCREEN
  end_message;
  
end

%% DEFINITION OF THE BODY FORCE FUNCTIONS FOR THE CANTILEVER BEAM

function bx = bx_body_force_fun(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the return "bx"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "bx" is an array that has the same form of the input "x" or "y"
  
  bx=0;  % is used as a force per volume 
end
function by = by_body_force_fun(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the return "by"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "by" is an array that has the same form of the input "x" or "y"
  
  by=0;  % is used as a force per volume  
end

%% DEFINITION OF NEUMANN FUNCTIONS FOR THE CANTILEVER BEAM

function fx = fx_Neumann_fun(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the return "fx"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "fx" is an array that has the same form of the input "x" or "y"
  
  fx=0;   % is used as a force per length 
end
function fy = fy_Neumann_fun(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the return "fy"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "fy" is an array that has the same form of the input "x" or "y"
  
  D=4; Ix=1*D*D*D/12; P=-1000;
  fy=P*(D*D/4 - y.*y)/(2*Ix);    % is used as a force per length 
end

%% DEFINITION OF DIRICHLET FUNCTIONS FOR THE CANTILEVER BEAM

function ux = ux_Dirichlet_fun(x,y,matProps)
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
  Ey=matProps.Ey; nu=matProps.nu; 
  if strcmp(matProps.plane_state,'plane_stress')
    Ey_bar=Ey; nu_bar=nu;
  elseif strcmp(matProps.plane_state,'plane_strain')
    Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  end
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  ux(:,1)=(-P/(6*Ey_bar*Ix))*y.*((6*L-3*x).*x +(2+nu_bar)*y.*y-3*D*D/2*(1+nu_bar));    
end
function uy = uy_Dirichlet_fun(x,y,matProps)
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
  Ey=matProps.Ey; nu=matProps.nu; 
  if strcmp(matProps.plane_state,'plane_stress')
    Ey_bar=Ey; nu_bar=nu;
  elseif strcmp(matProps.plane_state,'plane_strain')
    Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  end
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  uy(:,1)=P/(6*Ey_bar*Ix)*(3*nu_bar*y.*y.*(L-x)+(3*L-x).*x.*x);     
end

%% DEFINITION OF THE EXACT SOLUTIONS FOR THE CANTILEVER BEAM

function ux = ux_exact(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "ux" is an array that has the same form of the input "x" or "y"
  
  Ey=matProps.Ey; nu=matProps.nu; 
  if strcmp(matProps.plane_state,'plane_stress')
    Ey_bar=Ey; nu_bar=nu;
  elseif strcmp(matProps.plane_state,'plane_strain')
    Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  end
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  ux=(-P/(6*Ey_bar*Ix))*y.*((6*L-3*x).*x +(2+nu_bar)*y.*y-3*D*D/2*(1+nu_bar));    
end
function uy = uy_exact(x,y,matProps)  
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "uy" is an array that has the same form of the input "x" or "y"
  
  Ey=matProps.Ey; nu=matProps.nu; 
  if strcmp(matProps.plane_state,'plane_stress')
    Ey_bar=Ey; nu_bar=nu;
  elseif strcmp(matProps.plane_state,'plane_strain')
    Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  end
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  uy=P/(6*Ey_bar*Ix)*(3*nu_bar*y.*y.*(L-x)+(3*L-x).*x.*x);  
end
function p = p_exact(x,y,matProps)  
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "uy" is an array that has the same form of the input "x" or "y"
  
  P=-1000;
  D=4; 
  L=8; 
  Ix=1*D*D*D/12;
    
  fx=-P/Ix*(L-x).*y;
  fy=zeros(length(x),1);
  fxy=P/(2*Ix)*(D*D/4-y.^2);

  D=matProps.D; % defined as per Gain et al.
                % to come back to the std definiton in FEM books:
                % D(3,3)/4, this gives the following stress:
                % s12 = D(3,3)/4*(2*e12) = (D(3,3)/2)*e12
                % So, to obtain strain in the format [e11 e22 e12]
                % one needs to do: e12 = (D(3,3)/2)^-1 * s12
                
  D(3,3)=D(3,3)/2; % to obtain the strain as [e11 e22 e12]    
  strain_vec = D\[fx';fy';fxy'];
  
  p = -matProps.lam*(strain_vec(1,:)' + strain_vec(2,:)');   % update the value 
 
end

function strainvec = strainvec_exact(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "duxdx" is an array that has the same form of the input "x" or "y"
  
  Ey=matProps.Ey; nu=matProps.nu; 
  if strcmp(matProps.plane_state,'plane_stress')
    Ey_bar=Ey; nu_bar=nu;
  elseif strcmp(matProps.plane_state,'plane_strain')
    Ey_bar=Ey/(1-nu*nu); nu_bar=nu/(1-nu);
  end
  D=4; L=8; Ix=1*D*D*D/12; P=-1000;
  duxdx=-(P*y.*(6*L-6*x))/(6*Ey_bar*Ix);
  duydx=-(P*(3*nu_bar*y.^2-2*x.*(3*L-x)+x.^2))/(6*Ey_bar*Ix);  
  duxdy=-(P*((nu_bar+2)*y.^2+x.*(6*L-3*x)-(3*D^2*(nu_bar+1))/2))/(6*Ey_bar*Ix)...
        -(P*y.^2*(nu_bar+2))/(3*Ey_bar*Ix);  
  duydy=(P*nu_bar*y.*(L-x))/(Ey_bar*Ix);      
 
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

