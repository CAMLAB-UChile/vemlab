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
%                            L-shaped Domain
%                    H1 seminorm vs. degrees of fredom
%                           2D Poisson Problem
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

function L_shaped_domain_poisson2d_H1_vs_dof

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

  mesh_filename_polylla_random=...
                {'Lshape_poisson2d_104poly_elems_random_ssalinas.txt',...
                 'Lshape_poisson2d_1160poly_elems_random_ssalinas.txt',...
                 'Lshape_poisson2d_5523poly_elems_random_ssalinas.txt',...
                 'Lshape_poisson2d_53650poly_elems_random_ssalinas.txt'}; 
               
  mesh_filename_polylla_semiuniform=...
                {'Lshape_poisson2d_134poly_elems_semiuniform_ssalinas.txt',...
                 'Lshape_poisson2d_1050poly_elems_semiuniform_ssalinas.txt',...
                 'Lshape_poisson2d_5042poly_elems_semiuniform_ssalinas.txt',...
                 'Lshape_poisson2d_51494poly_elems_semiuniform_ssalinas.txt'};      
               
  mesh_filename_voronoi_random=...
                {'Lshape_poisson2d_300poly_elems_voronoi_random_ssalinas.txt',...
                 'Lshape_poisson2d_3590poly_elems_voronoi_random_ssalinas.txt',...
                 'Lshape_poisson2d_17000poly_elems_voronoi_random_ssalinas.txt',...
                 'Lshape_poisson2d_168000poly_elems_voronoi_random_ssalinas.txt'};                 
      
  mesh_filename_voronoi_semiuniform=...
                {'Lshape_poisson2d_224poly_elems_voronoi_semiuniform_ssalinas.txt',...
                 'Lshape_poisson2d_1607poly_elems_voronoi_semiuniform_ssalinas.txt',...
                 'Lshape_poisson2d_7513poly_elems_voronoi_semiuniform_ssalinas.txt',...
                 'Lshape_poisson2d_76218poly_elems_voronoi_semiuniform_ssalinas.txt'};  
  
  xcm=zeros(4,4);
  ycm=zeros(4,4);
  
  for i=1:length(mesh_filename_polylla_random)      
    [H1norm,cpu_time,ndofs]=...
      L_shaped_domain_poisson2d(opsystem,vemlab_root_dir,mesh_filename_polylla_random{i},'VEM2D');
    xcm(1,i)=sqrt(ndofs);
    ycm(1,i)=H1norm;    
  end

  for i=1:length(mesh_filename_polylla_semiuniform)      
    [H1norm,cpu_time,ndofs]=...
      L_shaped_domain_poisson2d(opsystem,vemlab_root_dir,mesh_filename_polylla_semiuniform{i},'VEM2D');
    xcm(2,i)=sqrt(ndofs);
    ycm(2,i)=H1norm;       
  end

  for i=1:length(mesh_filename_voronoi_random)      
    [H1norm,cpu_time,ndofs]=...
      L_shaped_domain_poisson2d(opsystem,vemlab_root_dir,mesh_filename_voronoi_random{i},'VEM2D');
    xcm(3,i)=sqrt(ndofs);
    ycm(3,i)=H1norm;       
  end

  for i=1:length(mesh_filename_voronoi_semiuniform)      
    [H1norm,cpu_time,ndofs]=...
      L_shaped_domain_poisson2d(opsystem,vemlab_root_dir,mesh_filename_voronoi_semiuniform{i},'VEM2D');
    xcm(4,i)=sqrt(ndofs);
    ycm(4,i)=H1norm;      
  end
   
  H1norm_plot(xcm,ycm); % plot H1 seminorm vs DOF
                                     
end

function H1norm_plot(xcm,ycm)

  figure;
  
  xc1=xcm(1,:);
  yc1=ycm(1,:);
  xc2=xcm(2,:);
  yc2=ycm(2,:);
  xc3=xcm(3,:);
  yc3=ycm(3,:);  
  xc4=xcm(4,:);
  yc4=ycm(4,:);  
  
  % Create plot
  loglog(xc1, yc1, '-d', 'Color', '[0.8500 0.3250 0.0980]', 'MarkerFaceColor', '[0.6500 0.1250 0.0680]', 'linewidth', 2.0, 'markersize', 6); hold on;
  loglog(xc2, yc2, '-o', 'Color', '[0 0.4470 0.7410]', 'MarkerFaceColor', '[0 0.2470 0.5410]', 'linewidth', 2.0, 'markersize', 6); hold on;
  loglog(xc3, yc3, '-s', 'Color', '[0.4660 0.6740 0.1880]', 'MarkerFaceColor', '[0.2660 0.4740 0.0880]', 'linewidth', 2.0, 'markersize', 6); hold on; 
  loglog(xc4, yc4, '-*', 'Color', '[0.4940 0.1840 0.5560]', 'MarkerFaceColor', '[0.4060 0.1240 0.4580]', 'linewidth', 2.0, 'markersize', 6); hold on;
  
  %
  % add a triangle with labels 1 and 1 for the rate : pick 4th point
  %

  % xc1=xc2=xc3=xc4
  x0 = 8.5*min(xc1); % the x-coord of the bottom left corner of the triangle (the starting point)
  y0 = mean([yc1,yc2,yc3,yc4])*0.07;   % the y-coord of the bottom left corner of the triangle (the starting point)
  
  scale = 300;  % to scale the triangle
  x_shift_bottom = 3.5;  % to center in the x-direction the number that goes on the bottom side of the triangle
  y_shift_bottom = 0.25; % to center in the y-direction the number that goes on the bottom side of the triangle
  x_shift_right = 0.18;  % to center in the x-direction the number that goes on the right side of the triangle
  y_shift_right = 6.0; % to center in the y-direction the number that goes on the right side of the triangle

  slope = 1;  % set the slope of the triangle

  ybar = y0 * exp( slope * log( (x0+scale)/x0 ));
  x = [x0 x0+scale x0 x0]; y = [y0 y0 ybar y0];
  plot(x,y,'k','linewidth', 1.1);
  text(x0+scale/x_shift_bottom,y0-y0*y_shift_bottom,'$$1$$','interpreter','latex','fontsize',12);
  text(x0-x0*x_shift_right,y0+ybar/y_shift_right,'$$1$$','interpreter','latex','fontsize',12);  

  %
  % format the axis
  %

  fsize = 12; % fontsize
  xlabel('(DOF)$$^{1/2}$$','interpreter','latex','fontsize', fsize); 
  ylabel('$$\|\mathbf{u}-\mathbf{u}_h\|_{H^1}/\|\mathbf{u}\|_{H^1}$$','interpreter','latex','fontsize', fsize);

  axis([10^1, 10^3, 10^-3, 10^0]);
  set(gca, 'fontsize', fsize+1);
  set(gca, 'xtick', 10.^(1:3));
  set(gca, 'ytick', 10.^(-3:0));
  set(gca,'TickLabelInterpreter','latex');

  ax = gca;
  ax.XMinorTick = 'on';
  ax.YMinorTick = 'on';   

  %
  % format the legend
  % 

  s1 = append('VEM',' ','(Polylla - random)');
  s2 = append('VEM',' ','(Polylla - semiuniform)');
  s3 = append('VEM',' ','(Voronoi - random)');
  s4 = append('VEM',' ','(Voronoi - semiuniform)');
  
  legend1 = legend(s1,s2,s3,s4,'location','northeast');    
  set(legend1,'FontSize',fsize-2);

end

function [H1rel,cpu_time,ndofs]=L_shaped_domain_poisson2d(opsystem,vemlab_root_dir,mesh_filename,vemlab_method)

  %% INPUT DATA 
  
  % set material parameters
  k=1;                               % conductivity (isotropic material)
  
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
  
  % module
  vemlab_module='Poisson';    % 'LinearElastostatics' or 'Poisson'
  
  % solver
  vemlab_solver='sparse';   % 'sparse' or 'dense'  
  
  %% PLOT AND OUTPUT OPTIONS
  % to setup plot and output options go to folder 'config' and 
  % set the parameters inside function 'plot_and_output_options.m'    
  
  %% VEMLAB CONFIGURATION
  
  stability_type = []; % not used---only one stability implemented for Poisson module
  config=config_vemlab(opsystem,vemlab_root_dir,mesh_filename,vemlab_module,...
                       vemlab_method,vemlab_solver,stability_type);
                     
  % deactivate matlab's plots if these were not deactivated in
  % function plot_and_ouput_options.m
  config.create_matlab_contour_plots='no';  
  config.create_matlab_exact_contour_plots='no';                       
                                         
  %% PRINT INIT MESSAGE TO SCREEN

  init_message(config);

  %% PREPROCESSING
  
  tic
  
  % read mesh
  domainMesh=read_mesh(config); 
  ndofs=size(domainMesh.coords,1); % total number of degrees of freedom
  
  % plot mesh  
%   plot_mesh2d(domainMesh,config);

  fprintf('\n');
  fprintf('Simulation started...\n');  
  
  %% CREATE POISSON MATERIAL 
  
  matProps=material_parameters_poisson(k);
  
  %% SOURCE TERM FUNCTIONS
  
  source_term_fun_values=@(x,y,matProps)source_term_fun(x,y,matProps);
  
  %% ASSEMBLY 
  
  fprintf('Assemblying element matrices...\n'); 
  [K_global,f_global]=assembly(domainMesh,config,matProps,source_term_fun_values);
  
  %% DIRICHLET BOUNDARY NODES/DOFS/FUNCTIONS
  % (see definition of functions at the bottom of this file)
  % The current problem has the entire boundary with Dirichlet BCs and
  % domain_type is 'Custom'. Therefore, we use
  % domainMesh.boundary_nodes.Dirichlet and domainMesh.boundary_dofs.Dirichlet
  
  Dirichet_boundary_nodes=domainMesh.boundary_nodes.Dirichlet;
  Dirichet_boundary_dofs=domainMesh.boundary_dofs.Dirichlet; 
  Dirichlet_fun_values=@(x,y,matProps)u_Dirichlet_fun(x,y,matProps);   

  %% ENFORCE DIRICHLET BCS ON THE BOUNDARY NODES
  
  fprintf('Enforcing Dirichlet boundary conditions...\n');   
  u_nodal_sol=zeros(length(domainMesh.coords),1);    
  DB_dofs=compute_Dirichlet_BCs(domainMesh,config,Dirichet_boundary_nodes,...
                                Dirichet_boundary_dofs,Dirichlet_fun_values,matProps); % DOFs with Dirichlet BCs
  u_nodal_sol(DB_dofs.indexes)=DB_dofs.values;
  f_global=f_global-K_global*u_nodal_sol;    
  
  %% SOLVE RESULTING SYSTEM OF LINEAR EQUATIONS 
  
  fprintf('Solving system of linear equations...\n');   
  num_nodes=length(domainMesh.coords(:,1));    
  free_dofs=setdiff(1:num_nodes,DB_dofs.indexes); % dofs that do not have Dirichlet boundary conditions associated  
  u_nodal_sol(free_dofs)=K_global(free_dofs,free_dofs)\f_global(free_dofs); % nodal solution at free dofs
  
  cpu_time = toc;
  
  %% EXACT SOLUTIONS
  % (see definition of functions at the bottom of this file)
  
  exact_solution_handle.u=@(x,y,matProps)u_exact(x,y,matProps);
  exact_solution_handle.dudx=@(x,y,matProps)dudx_exact(x,y,matProps);   
  exact_solution_handle.dudy=@(x,y,matProps)dudy_exact(x,y,matProps);    
  
  %% NORMS OF THE ERROR
  
  fprintf('\n');
  fprintf('Computing norms of the solution error...\n');   
  [h_max,L2rel,H1rel,~]=compute_norms_of_the_error(exact_solution_handle,u_nodal_sol,domainMesh,config,matProps);  
  
%   %% POSTPROCESSING
%   
%   % plot numerical solutions  
%   postprocess_numerical_solution_poisson2d(domainMesh,u_nodal_sol,matProps,config);
%   
%   % plot exact solutions    
%   postprocess_exact_solution_poisson2d(domainMesh,exact_solution_handle,matProps,config);  
  
  %% PRINT END MESSAGE TO SCREEN
  end_message;  
  
end

%% DEFINITION OF THE SOURCE TERM FUNCTION

function b = source_term_fun(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the return "b"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "b" is an array that has the same form of the input "x" or "y"
  
  b=0.0;  % is used as a force per volume
end

%% DEFINITION OF DIRICHLET FUNCTIONS FOR THE SQUARE PLATE WITH HEAT SOURCE

function u = u_Dirichlet_fun(x,y,matProps)
  % INPUT: x,y are vectors containing the coordinates of the nodes lying on the 
  % Dirichlet boundary, therefore if the Dirichlet conditions depend on x and y,
  % consider using something like x.*y (i.e., use the dot symbol).
  %
  % OUTPUT: u = array in which its first column contains the value of the 
  % Dirichlet boundary condition for the degrees of freedom, and its second 
  % column = free (1) or fixed (0) to indicate whether the corresponding value 
  % of the degree of freedom in the first column should be ignored (free) or 
  % applied (fixed).
  
  N=length(x);
  u=zeros(N,2);  % first column = value; second column = free (1) or fixed (0)
                 % In this case, all dofs are fixed and have the values of u(:,1)
  r=sqrt(x.*x + y.*y);
  theta=atan2(y,x);
  % update the values
  ind=find(theta<0); % since in the interval (pi,2*pi) atan2(y,x) returns values in the interval (-pi,0), convert to interval (-pi,0) to (pi,2*pi)
  theta(ind)=pi()-(-pi()-theta(ind));
  u(:,1)=(r.^(2/3)).*sin(2/3*theta); 
end

%% DEFINITION OF THE EXACT SOLUTIONS

function u = u_exact(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "u" is an array that has the same form of the input "x" or "y"  
  
  r=sqrt(x.*x + y.*y);
  theta=atan2(y,x);
  ind=find(theta<0); % since in the interval (pi,2*pi) atan2(y,x) returns values in the interval (-pi,0), convert to interval (-pi,0) to (pi,2*pi)
  theta(ind)=pi()-(-pi()-theta(ind));
  u(:,1)=(r.^(2/3)).*sin(2/3*theta);
  
end

function dudx = dudx_exact(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "dudx" is an array that has the same form of the input "x" or "y" 
  
  theta=atan2(y,x);  
  ind=find(theta<0); % since in the interval (pi,2*pi) atan2(y,x) returns values in the interval (-pi,0), convert to interval (-pi,0) to (pi,2*pi)
  theta(ind)=pi()-(-pi()-theta(ind));  
  dudx=(2*(x.*sin(2/3*theta) - y.*cos(2/3*theta)))./(3*(x.*x + y.*y).^(2/3));
end

function dudy = dudy_exact(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "dudy" is an array that has the same form of the input "x" or "y"  
  
  theta=atan2(y,x);  
  ind=find(theta<0); % since in the interval (pi,2*pi) atan2(y,x) returns values in the interval (-pi,0), convert to interval (-pi,0) to (pi,2*pi)
  theta(ind)=pi()-(-pi()-theta(ind)); 
  dudy=(2*(y.*sin(2/3*theta) + x.*cos(2/3*theta)))./(3*(x.*x + y.*y).^(2/3));
end
