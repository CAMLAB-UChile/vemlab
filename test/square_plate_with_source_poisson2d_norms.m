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
%                   Square Plate (1x1) with Heat Source
%                        L^2 norm and H^1 seminorm
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

function square_plate_with_source_poisson2d_norms

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
               
  mesh_filename={'square_plate_poisson2d_20poly_elems.txt',...
                 'square_plate_poisson2d_100poly_elems.txt',...
                 'square_plate_poisson2d_1000poly_elems.txt',...
                 'square_plate_poisson2d_6000poly_elems.txt'};                  
                                      
  xcm=zeros(1,4);
  ycm=zeros(1,4);
  ycm1=zeros(1,4);    
  for i=1:length(mesh_filename)      
    [h_max_1,L2rel_1,H1rel_1]=...
      square_plate_with_source_poisson2d(opsystem,vemlab_root_dir,mesh_filename{i},'VEM2D');
    xcm(1,i)=h_max_1;
    ycm(1,i)=H1rel_1;
    ycm1(1,i)=L2rel_1;      
  end
  textnote='';% 'some text', '$\lambda/\mu=5\times 10^7$', etc.
  H1plot(xcm,ycm,textnote); % plot H1 norm
  L2plot(xcm,ycm1,textnote); % plot L2 norm 
                                     
end

function H1plot(xcm,ycm,textnote)

  figure;
  
  xc=xcm(1,:);
  yc1=ycm(1,:);

  % Create plot
  
  loglog(xc, yc1, '-d', 'Color', '[0.8500 0.3250 0.0980]', 'MarkerFaceColor', '[0.6500 0.1250 0.0680]', 'linewidth', 2.0, 'markersize', 6); hold on;

  %
  % add a triangle with labels 1 and 1 for the rate : pick 4th point
  %

  x0 = min(xc); % the x-coord of the bottom left corner of the triangle (the starting point) 
  y0 = min(yc1)*0.5;   % the y-coord of the bottom left corner of the triangle (the starting point) 

  scale = 0.015;  % to scale the triangle
  x_shift_bottom = 3.0;  % to center in the x-direction the number that goes on the bottom side of the triangle
  y_shift_bottom = 0.3; % to center in the y-direction the number that goes on the bottom side of the triangle
  x_shift_right = 0.06;  % to center in the x-direction the number that goes on the right side of the triangle
  y_shift_right = 5.0; % to center in the y-direction the number that goes on the right side of the triangle

  slope = 1;  % set the slope of the triangle

  ybar = y0 * exp( slope * log( (x0+scale)/x0 )) ;
  x = [x0 x0+scale x0+scale x0]; y = [y0 y0 ybar y0];
  plot(x,y,'k','linewidth', 1.6);
  text(x0+scale/x_shift_bottom,y0-y0*y_shift_bottom,'$$1$$','interpreter','latex','fontsize',12);
  text(x0+scale+(x0+scale)*x_shift_right,y0+ybar/y_shift_right,'$$1$$','interpreter','latex','fontsize',12);

  %
  % format the axis
  %

  fsize = 12; % fontsize
  xlabel('Element size','interpreter','latex','fontsize', fsize); 
  ylabel('$$\|\mathbf{u}-\mathbf{u}_h\|_{H^1}/\|\mathbf{u}\|_{H^1}$$','interpreter','latex','fontsize', fsize);

  axis([10^-3, 10^0, 10^-3, 10^0]);
  set(gca, 'fontsize', fsize+1);
  set(gca, 'xtick', 10.^(-3:0));
  set(gca, 'ytick', 10.^(-3:0));
  set(gca,'TickLabelInterpreter','latex');

  ax = gca;
  ax.XMinorTick = 'on';
  ax.YMinorTick = 'on';   

  %
  % format the legend
  % 
  
  s1 = 'VEM';
%   s1 = append('VEM',' ','(',textnote,')');

  legend1 = legend(s1,'location','southeast');    
  %set(legend1,'FontSize',fsize,'interpreter','latex');  
  set(legend1,'FontSize',fsize-2);
  
%   annotation('textbox',...
%       [0.154571428571427 0.830952382158672 0.139285710666861 0.0642857130794299],...
%       'String',textnote,'FitBoxToText','on','Interpreter','latex');   

end

function L2plot(xcm,ycm,textnote)

  figure;
  
  xc=xcm(1,:);
  yc1=ycm(1,:);

  % Create plot
  
  loglog(xc, yc1, '-d', 'Color', '[0.8500 0.3250 0.0980]', 'MarkerFaceColor', '[0.6500 0.1250 0.0680]', 'linewidth', 2.0, 'markersize', 6); hold on;

  %
  % add a triangle with labels 1 and 1 for the rate : pick 4th point
  %

  x0 = min(xc); % the x-coord of the bottom left corner of the triangle (the starting point)
  y0 = min(yc1)*0.25;   % the y-coord of the bottom left corner of the triangle (the starting point)    

  scale = 0.015;  % to scale the triangle
  x_shift_bottom = 3.0;  % to center in the x-direction the number that goes on the bottom side of the triangle
  y_shift_bottom = 0.48; % to center in the y-direction the number that goes on the bottom side of the triangle
  x_shift_right = 0.06;  % to center in the x-direction the number that goes on the right side of the triangle
  y_shift_right = 7.0; % to center in the y-direction the number that goes on the right side of the triangle

  slope = 2;  % set the slope of the triangle

  ybar = y0 * exp( slope * log( (x0+scale)/x0 )) ;
  x = [x0 x0+scale x0+scale x0]; y = [y0 y0 ybar y0];
  plot(x,y,'k','linewidth', 1.6);
  text(x0+scale/x_shift_bottom,y0-y0*y_shift_bottom,'$$1$$','interpreter','latex','fontsize',12);
  text(x0+scale+(x0+scale)*x_shift_right,y0+ybar/y_shift_right,'$$2$$','interpreter','latex','fontsize',12);

  %
  % format the axis
  %

  fsize = 12; % fontsize
  xlabel('Element size','interpreter','latex','fontsize', fsize); 
  ylabel('$$\|\mathbf{u}-\mathbf{u}_h\|_{L^2}/\|\mathbf{u}\|_{L^2}$$','interpreter','latex','fontsize', fsize);

  axis([10^-3, 10^0, 10^-5, 10^-1]);
  set(gca, 'fontsize', fsize+1);
  set(gca, 'xtick', 10.^(-3:0));
  set(gca, 'ytick', 10.^(-5:-1));
  set(gca,'TickLabelInterpreter','latex');

  ax = gca;
  ax.XMinorTick = 'on';
  ax.YMinorTick = 'on';   

  %
  % format the legend
  % 
  
  s1 = 'VEM';
%   s1 = append('VEM',' ','(',textnote,')');
  
  legend1 = legend(s1,'location','southeast');   
  %set(legend1,'FontSize',fsize,'interpreter','latex');
  set(legend1,'FontSize',fsize-2);
  
%   annotation('textbox',...
%       [0.154571428571427 0.830952382158672 0.139285710666861 0.0642857130794299],...
%       'String',textnote,'FitBoxToText','on','Interpreter','latex'); 

end

function [h_max,L2rel,H1rel]=square_plate_with_source_poisson2d(opsystem,vemlab_root_dir,mesh_filename,vemlab_method)

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
  
  % read mesh
  domainMesh=read_mesh(config);
  
  % plot mesh  
  plot_mesh2d(domainMesh,config);

  fprintf('\n');
  fprintf('Simulation started...\n');  
  
  %% CREATE POISSON MATERIAL 
  
  matProps=material_parameters_poisson(k);
  
  %% SOURCE TERM FUNCTIONS
  
  source_term_fun_values=@(x,y,matProps)source_term_fun(x,y,matProps);
  
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
  toc
  
  %% EXACT SOLUTIONS
  % (see definition of functions at the bottom of this file)
  
  exact_solution_handle.u=@(x,y,matProps)u_exact(x,y,matProps);
  exact_solution_handle.dudx=@(x,y,matProps)dudx_exact(x,y,matProps);   
  exact_solution_handle.dudy=@(x,y,matProps)dudy_exact(x,y,matProps);    
  
  %% NORMS OF THE ERROR
  
  fprintf('\n');
  fprintf('Computing norms of the solution error...\n');   
  [h_max,L2rel,H1rel,~]=compute_norms_of_the_error(exact_solution_handle,u_nodal_sol,domainMesh,config,matProps);  
  
  %% POSTPROCESSING
  
  % plot numerical solutions  
  postprocess_numerical_solution_poisson2d(domainMesh,u_nodal_sol,matProps,config);
  
  % plot exact solutions    
  postprocess_exact_solution_poisson2d(domainMesh,exact_solution_handle,matProps,config);  
  
  %% PRINT END MESSAGE TO SCREEN
  end_message;  
  
end

%% DEFINITION OF THE SOURCE TERM FUNCTION

function b = source_term_fun(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the return "b"
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "b" is an array that has the same form of the input "x" or "y"
  
  b=32*y.*(1-y) + 32*x.*(1-x);  % is used as a force per volume
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
end

%% DEFINITION OF THE EXACT SOLUTIONS

function u = u_exact(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "u" is an array that has the same form of the input "x" or "y"  
  
  u=16*x.*y.*(1-x).*(1-y);
end

function dudx = dudx_exact(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "dudx" is an array that has the same form of the input "x" or "y"  
  
  dudx=16*y.*(1-y).*(1-2*x);
end

function dudy = dudy_exact(x,y,matProps)
  % Use something like x.*y (i.e., use the dot symbol) if the exact solution
  % depends on x and y. This way, this function will also serve in case x and y 
  % are arrays. If the function does not depend on x and y, make sure that the
  % return "dudy" is an array that has the same form of the input "x" or "y"  
  
  dudy=16*x.*(1-x).*(1-2*y);
end