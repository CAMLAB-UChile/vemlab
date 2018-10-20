function [create_matlab_contour_plots,plot_mesh,plot_mesh_over_results,write_solutions_to_text_file,...
          write_solutions_to_GiD_file,write_solutions_to_VTK_file,poisson2d_plot_scalar_field,...
          poisson2d_plot_flux,poisson2d_plot_grad,linelast2d_plot_displacement,...
          linelast2d_plot_stress,linelast2d_plot_strain,linelast2d_plot_deformed_domain,...
          linelast2d_scale_for_plotting_deformed_domain] = plot_and_output_options
  
  % INSTRUCTIONS: SET VARIABLES WITH 'yes' or 'no' ONLY, except 
  %               linelast2d_scale_for_plotting_deformed_domain, which must be
  %               assigned a positive number.
  %
  % NOTE: Currently, the MATLAB plotting of stresses and strains 
  %       (linear elastostatic problem) or fluxes and gradients (Poisson problem)
  %       is very limited and won't recognize holes that might come with
  %       the domain geometries. It will work ok with rectangular geometries
  %       without holes, but the mesh will need to be very refined to obtain a
  %       colormap on the whole geometry ... just try a coarse mesh and then a
  %       refined mesh to see what is being said.
  %
  %       If quality colormap plots are required for stresses, strains, fluxes
  %       and gradients, it is highly recommended writing the GiD output files and
  %       visualizing them in "GiD the pre and postprocessor" (www.gidhome.com)
  %
  disp('*** DONT FORGET TO SETUP THE PLOT AND OUTPUT OPTIONS ***');
  disp(' ');
  disp('To setup them go to file config\plot_and_output_options.m');
  disp(' ');  
  disp('NOTE: Currently, the MATLAB plotting of stresses and strains'); 
  disp('(linear elastostatic problem) or fluxes and gradients (Poisson problem)');
  disp('is very limited and wont recognize holes that might come with');
  disp('the domain geometries. It will work fine with rectangular geometries');
  disp('without holes, but the mesh will need to be very refined to obtain');
  disp('quality colormaps on the whole geometry ... just try a coarse mesh');
  disp('and then a refined mesh to see what is being said.');
  disp(' ');    
  disp('If quality colormap plots are required for stresses, strains, fluxes');
  disp('and gradients, it is higly recommended writing the GiD output files and');
  disp('visualizing them in GiD the pre and postprocessor (www.gidhome.com)');
  disp(' '); 
  disp('CAUTION: however, a correct visualization in GiD requires a poygonal mesh');
  disp('that is formed by CONVEX POLYGONS.');
  disp(' ');   
  disp('Press any key to continue ...');  
  disp(' ');    
  pause;
  
  %% GENERAL
  
  create_matlab_contour_plots='yes';
  plot_mesh='yes';
  plot_mesh_over_results='yes';
  write_solutions_to_text_file='yes';
  write_solutions_to_GiD_file='yes';
  write_solutions_to_VTK_file='yes';  
  
  %% POISSON MODULE
  
  % plotting of main variables to MATLAB/GiD/VTK figures
  poisson2d_plot_scalar_field.u='yes';
  
  % plotting of fluxes to MATLAB/GiD figures
  poisson2d_plot_flux.qx='no';
  poisson2d_plot_flux.qy='no'; 
  poisson2d_plot_flux.qnorm='yes';   % norm of the flux
  
  % plotting of gradients to MATLAB/GiD figures
  poisson2d_plot_grad.dx='no';
  poisson2d_plot_grad.dy='no';   
  poisson2d_plot_grad.dnorm='yes';   % norm of the gradient
  
  %% LINELAST MODULE
  
  % options for plotting deformed domain to MATLAB figures
  linelast2d_plot_deformed_domain='yes';
  linelast2d_scale_for_plotting_deformed_domain=1; % a number > 1 will scale the deformed domain when plotting to MATLAB figures

  % plotting of main variables to MATLAB/GiD/VTK figures
  linelast2d_plot_displacement.ux='yes';  
  linelast2d_plot_displacement.uy='yes';  
  linelast2d_plot_displacement.unorm='yes';  % norm of the displacement
  
  % plotting of stresses to MATLAB/GiD figures
  linelast2d_plot_stress.s11='no';
  linelast2d_plot_stress.s12='no';
  linelast2d_plot_stress.s22='no';
  linelast2d_plot_stress.s33='no';
  linelast2d_plot_stress.s1='no';
  linelast2d_plot_stress.s2='no';  
  linelast2d_plot_stress.s3='no';  
  linelast2d_plot_stress.vm='no';  
  
  % plotting of strains to MATLAB/GiD figures
  linelast2d_plot_strain.e11='no';
  linelast2d_plot_strain.e12='no';
  linelast2d_plot_strain.e22='no';
  linelast2d_plot_strain.e33='no';
  linelast2d_plot_strain.e1='no';
  linelast2d_plot_strain.e2='no';  
  linelast2d_plot_strain.e3='no';   
  
end

