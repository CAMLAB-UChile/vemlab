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
%                        L^2 norm and H^1 seminorm
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

function cantilever_beam_linelast2d_norms

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
  
  mesh_filename={'cantilever_beam_250poly_elems.txt',...
                 'cantilever_beam_1500poly_elems.txt',...
                 'cantilever_beam_12000poly_elems.txt',...
                 'cantilever_beam_96000poly_elems.txt'};          
               
  % Stability
  stability_type = [1,2,3];     % 0 -> without stability
                                % 1 -> Gain et al.
                                % 2 -> D-recipe
                                % 3 -> modified D-recipe                     
                              
  % E=1e7, nu=0.3 
  Ey=10^7;                               % Young's modulus
  nu=0.3;                                % Poisson's ratio
  elastic_state='plane_strain';          % 'plane_strain' or 'plane_stress'    
  xcm=zeros(4,4);
  ycm=zeros(4,4);
  ycm1=zeros(4,4);  
  ycm2=zeros(4,4);   
  for i=1:length(mesh_filename)
    for j=1:length(stability_type)                            
      [h_max_1,L2rel_1,H1rel_1,L2prel_1]=...
        cantilever_beam_linelast2d(opsystem,vemlab_root_dir,mesh_filename{i},stability_type(j),'VEM2D',Ey,nu,elastic_state);
      xcm(j,i)=h_max_1;
      ycm(j,i)=H1rel_1;
      ycm1(j,i)=L2rel_1;
      ycm2(j,i)=L2prel_1;      
    end   
  end
  textnote = '$\nu = 0.3$';
  H1plot(xcm,ycm,textnote); % plot H1 norm
  L2plot(xcm,ycm1,textnote); % plot L2 norm
  L2pressureplot(xcm,ycm2,textnote); % plot L2 pressure norm
                                     
end

function H1plot(xcm,ycm,textnote)

  figure;
  
  xc=xcm(1,:);
  yc1=ycm(1,:);
  yc3=ycm(2,:);
  yc5=ycm(3,:);

  % Create plot
  loglog(xc, yc1, '-d', 'Color', '[0.8500 0.3250 0.0980]', 'MarkerFaceColor', '[0.6500 0.1250 0.0680]', 'linewidth', 2.0, 'markersize', 6); hold on;
  loglog(xc, yc3, '-o', 'Color', '[0 0.4470 0.7410]', 'MarkerFaceColor', '[0 0.2470 0.5410]', 'linewidth', 2.0, 'markersize', 6); hold on;
  loglog(xc, yc5, '-s', 'Color', '[0.4660 0.6740 0.1880]', 'MarkerFaceColor', '[0.2660 0.4740 0.0880]', 'linewidth', 2.0, 'markersize', 6); hold on;

  %
  % add a triangle with labels 1 and 1 for the rate : pick 4th point
  %

  x0 = min(xc); % the x-coord of the bottom left corner of the triangle (the starting point)
  y0 = min([yc1,yc3,yc5])*0.5;   % the y-coord of the bottom left corner of the triangle (the starting point)  

  scale = 0.1;  % to scale the triangle
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

  axis([10^-2, 10^0, 10^-4, 10^0]);
  set(gca, 'fontsize', fsize+1);
  set(gca, 'xtick', 10.^(-2:0));
  set(gca, 'ytick', 10.^(-4:0));
  set(gca,'TickLabelInterpreter','latex');

  ax = gca;
  ax.XMinorTick = 'on';
  ax.YMinorTick = 'on'; 

  %
  % format the legend
  % 

  s1 = 'VEM, Stab. 1';
  s3 = 'VEM, Stab. 2';
  s5 = 'VEM, Stab. 3';

  %legend1 = legend(s1,s3,s5,'Interpreter','latex','location','southeast');
  legend1 = legend(s1,s3,s5,'Interpreter','latex','Position',[0.659060736582002 0.146493338064267 0.227248787227522 0.184459042888114]);
  set(legend1,'FontSize',fsize-2);
  
  annotation('textbox',...
      [0.154571428571427 0.830952382158672 0.139285710666861 0.0642857130794299],...
      'String',textnote,'FitBoxToText','on','Interpreter','latex'); 

end

function L2plot(xcm,ycm,textnote)

  figure;
  
  xc=xcm(1,:);
  yc1=ycm(1,:);
  yc3=ycm(2,:);
  yc5=ycm(3,:);

  % Create plot
  loglog(xc, yc1, '-d', 'Color', '[0.8500 0.3250 0.0980]', 'MarkerFaceColor', '[0.6500 0.1250 0.0680]', 'linewidth', 2.0, 'markersize', 6); hold on;
  loglog(xc, yc3, '-o', 'Color', '[0 0.4470 0.7410]', 'MarkerFaceColor', '[0 0.2470 0.5410]', 'linewidth', 2.0, 'markersize', 6); hold on;
  loglog(xc, yc5, '-s', 'Color', '[0.4660 0.6740 0.1880]', 'MarkerFaceColor', '[0.2660 0.4740 0.0880]', 'linewidth', 2.0, 'markersize', 6); hold on;

  %
  % add a triangle with labels 1 and 1 for the rate : pick 4th point
  %

  x0 = min(xc); % the x-coord of the bottom left corner of the triangle (the starting point)
  y0 = min([yc1,yc3,yc5])*0.25;   % the y-coord of the bottom left corner of the triangle (the starting point)  

  scale = 0.1;  % to scale the triangle
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

  axis([10^-2, 10^0, 10^-7, 10^0]);
  set(gca, 'fontsize', fsize+1);
  set(gca, 'xtick', 10.^(-2:0));
  set(gca, 'ytick', 10.^(-7:0));
  set(gca,'TickLabelInterpreter','latex');

  ax = gca;
  ax.XMinorTick = 'on';
  ax.YMinorTick = 'on';   

  %
  % format the legend
  % 

  s1 = 'VEM, Stab. 1';
  s3 = 'VEM, Stab. 2';
  s5 = 'VEM, Stab. 3';

  %legend1 = legend(s1,s3,s5,'Interpreter','latex','location','southeast');
  legend1 = legend(s1,s3,s5,'Interpreter','latex','Position',[0.659060736582002 0.146493338064267 0.227248787227522 0.184459042888114]);
  set(legend1,'FontSize',fsize-2);
  
  annotation('textbox',...
      [0.154571428571427 0.830952382158672 0.139285710666861 0.0642857130794299],...
      'String',textnote,'FitBoxToText','on','Interpreter','latex'); 

end

function L2pressureplot(xcm,ycm,textnote)

  figure;
  
  xc=xcm(1,:);
  yc1=ycm(1,:);
  yc3=ycm(2,:);
  yc5=ycm(3,:);

  % Create plot
  loglog(xc, yc1, '-d', 'Color', '[0.8500 0.3250 0.0980]', 'MarkerFaceColor', '[0.6500 0.1250 0.0680]', 'linewidth', 2.0, 'markersize', 6); hold on;
  loglog(xc, yc3, '-o', 'Color', '[0 0.4470 0.7410]', 'MarkerFaceColor', '[0 0.2470 0.5410]', 'linewidth', 2.0, 'markersize', 6); hold on;
  loglog(xc, yc5, '-s', 'Color', '[0.4660 0.6740 0.1880]', 'MarkerFaceColor', '[0.2660 0.4740 0.0880]', 'linewidth', 2.0, 'markersize', 6); hold on;

  %
  % add a triangle with labels 1 and 1 for the rate : pick 4th point
  %

  x0 = min(xc); % the x-coord of the bottom left corner of the triangle (the starting point)
  y0 = min([yc1,yc3,yc5])*0.5;   % the y-coord of the bottom left corner of the triangle (the starting point)  

  scale = 0.1;  % to scale the triangle
  x_shift_bottom = 3.0;  % to center in the x-direction the number that goes on the bottom side of the triangle
  y_shift_bottom = 0.3; % to center in the y-direction the number that goes on the bottom side of the triangle
  x_shift_right = 0.06;  % to center in the x-direction the number that goes on the right side of the triangle
  y_shift_right = 5.0; % to center in the y-direction the number that goes on the right side of the triangle

  slope = 1.0;  % set the slope of the triangle

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
  ylabel('$$\|p-p_h\|_{L^2}/\|p\|_{L^2}$$','interpreter','latex','fontsize', fsize);

  axis([10^-2, 10^0, 10^-4, 10^1]);
  set(gca, 'fontsize', fsize+1);
  set(gca, 'xtick', 10.^(-2:0));
  set(gca, 'ytick', 10.^(-4:1));
  set(gca,'TickLabelInterpreter','latex');

  ax = gca;
  ax.XMinorTick = 'on';
  ax.YMinorTick = 'on';   

  %
  % format the legend
  % 

  s1 = 'VEM, Stab. 1';
  s3 = 'VEM, Stab. 2';
  s5 = 'VEM, Stab. 3';

  %legend1 = legend(s1,s3,s5,'Interpreter','latex','location','southeast');
  legend1 = legend(s1,s3,s5,'Interpreter','latex','Position',[0.659060736582002 0.146493338064267 0.227248787227522 0.184459042888114]);
  set(legend1,'FontSize',fsize-2);
  
  annotation('textbox',...
      [0.154571428571427 0.830952382158672 0.139285710666861 0.0642857130794299],...
      'String',textnote,'FitBoxToText','on','Interpreter','latex'); 

end

function [h_max,L2rel,H1rel,L2prel]=cantilever_beam_linelast2d(opsystem,vemlab_root_dir,mesh_filename,stability_type,vemlab_method,Ey,nu,elastic_state)

  %% INPUT DATA 
  
  % mesh filename: see available sample mesh files in folder "test/mesh_files/"
  %
  %     VEM2D meshes have the keyword "poly" in the mesh filename
  %     FEM2DQ4 meshes have the keyword "q4" in the mesh filename
  %     FEM2DT3 meshes have the keyword "t3" in the mesh filename 
  % 
  %     VEM2D meshes can be used with 'VEM2D' method only.
  %     FEM2DQ4 meshes can be used with 'FEM2DQ4' and 'VEM2D' methods.
  %     FEM2DT3 meshes can be used with 'FEM2DT3' and 'VEM2D' methods.
  

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
  [h_max,L2rel,H1rel,L2prel]=compute_norms_of_the_error(exact_solution_handle,u_nodal_sol,domainMesh,config,matProps);  
  
  %% POSTPROCESSING
  
  % plot numerical solutions  
  %postprocess_numerical_solution_linelast2d(domainMesh,u_nodal_sol,matProps,config);
  
  % plot exact solutions    
  %postprocess_exact_solution_linelast2d(domainMesh,exact_solution_handle,matProps,config); 
  
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
  
%   p
%   p2 = P/(2*Ix)*(L-x).*y
%   p3 = -(fx+fy+fz)/3

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