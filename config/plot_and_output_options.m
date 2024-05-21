function poo = plot_and_output_options
  
  % INSTRUCTIONS: SET VARIABLES WITH 'yes' or 'no', or numbers where
  % appropriate
  %
  disp('*** DONT FORGET TO SETUP THE PLOT AND OUTPUT OPTIONS ***');
  disp(' ');
  disp('To setup them go to file config\plot_and_output_options.m');
  disp(' ');  
  disp('Execution will resume in 4 seconds...');  
  disp(' ');    
  pause(4);
  
  %% GENERAL
  
  %---------------------------------------------------------
  % 1. Delete existing output files
  %---------------------------------------------------------
  
  poo.delete_existing_output_figures='yes'; % 'yes' or 'no'
  poo.delete_existing_GiD_output_files='yes'; % 'yes' or 'no'
  poo.delete_existing_txt_output_files='yes'; % 'yes' or 'no'
  poo.delete_existing_VTK_output_files='yes'; % 'yes' or 'no'
    
  %---------------------------------------------------------
  % 2. Plot and figures settings
  %---------------------------------------------------------
  
  % 2.1. Mesh
  poo.plot_mesh='yes';  
  poo.plot_mesh_linewidth = 1.0;  
  poo.plot_mesh_nodes = 'no';
  poo.plot_mesh_nodesize = 5.0;
  poo.plot_mesh_axis = 'yes';    
  
  % 2.2. General settings for figures
  poo.plot_vemlab_logo='yes'; % 'yes' or 'no'
  poo.create_matlab_contour_plots='yes'; % 'yes' or 'no' (for numerical solution)
  poo.create_matlab_exact_contour_plots='yes'; % 'yes' or 'no' (for exact solution)
  poo.front_matlab_plot_view_orientation='yes'; % 'yes' for front orientation, otherwise the plot view is isometric
  poo.plot_mesh_over_results='yes'; % 'yes' or 'no'
  poo.plot_figure_axis='yes'; % 'yes' or 'no'
  poo.print_figures='yes'; % 'yes' or 'no' (for numerical solution)
  poo.print_exact_figures='yes'; % 'yes' or 'no' (for exact solution)
  poo.save_matlab_figures='no'; % 'yes' or 'no' (for numerical solution) 
  poo.save_exact_matlab_figures='no'; % 'yes' or 'no' (for exact solution)   
  poo.print_figures_resolution=600; % 96, 150, 300, 600, 1200, 2400, etc.... any real number > 0
  
  % 2.3. Color system for the colorbar
  poo.matlab_colormap='RdYlBu';   % 'RdYlBu', 'spectral', 'parula', 'jet', 'hsv', 'hot', 'cool', 'spring', 'summer', 'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', 'default'
  poo.colormap_number_of_colors=[]; % leave an empty [] for matlab's default value, or specify an integer value between 1 and 256.
  
  % 2.4. Colorbar ticks and lines (box and ticks)
  poo.colorbar_tick_label_notation='%.2e';    % '%.2e', '%.2f', '%.4e', '%.4f', etc
  poo.colorbar_TickLength=0.013; % 0.01, 0.02, 0.03, etc.
  poo.colorbar_LineWidth=1.05; % 0.5, 1.0, 1.1, etc.
  
  % 2.5. Figure ticks
  poo.axis_xtick=[]; % leave [] for defaults xtick or specify like [0 2 4 6] or [0.1 0.4 0.8 1.2 1.6 2.0], etc.
  poo.axis_ytick=[]; % leave [] for defaults xtick or specify like [0 2 4 6] or [0.1 0.4 0.8 1.2 1.6 2.0], etc.
  poo.axis_ztick=[]; % leave [] for defaults xtick or specify like [0 2 4 6] or [0.1 0.4 0.8 1.2 1.6 2.0], etc. 
  
  poo.figure_XMinorTick = 'on'; % 'on' or 'off'
  poo.figure_YMinorTick = 'on'; % 'on' or 'off'
  poo.figure_ZMinorTick = 'on'; % 'on' or 'off'
  
  % 2.6. Figure grid
  poo.figure_XGrid = 'off'; % 'on' or 'off'
  poo.figure_YGrid = 'off'; % 'on' or 'off'
  poo.figure_ZGrid = 'off'; % 'on' or 'off'
  poo.figure_XMinorGrid = 'off'; % 'on' or 'off'
  poo.figure_YMinorGrid = 'off'; % 'on' or 'off'    
  poo.figure_ZMinorGrid = 'off'; % 'on' or 'off'  
  
  % 2.7. Setting grid and minor grid color to 'none' results in grid off
  poo.figure_GridColor = 'default'; % 'default', 'red', 'green', blue', 'cyan', 'magenta', 'yellow', 'black', 'white', 'none'
  poo.figure_MinorGridColor = 'default'; % 'default', 'red', 'green', blue', 'cyan', 'magenta', 'yellow', 'black', 'white', 'none'
  
  %---------------------------------------------------------  
  % 3. result output files
  %---------------------------------------------------------  
  
  poo.write_solutions_to_text_file='no';
  poo.write_solutions_to_GiD_file='no';
  poo.write_solutions_to_VTK_file='no';  
  poo.write_solutions_to_CPP_file='no';  
  
  
  %% POISSON MODULE
  
  % plotting of main variables to MATLAB/GiD/VTK figures
  poo.poisson2d_plot_scalar_field.u='yes';
  % colorbar limits  
  poo.poisson2d_plot_scalar_field.clim.u=[]; % leave an empty [] for matlab's default value, or specify [min max] values.   
  
  % plotting of fluxes to MATLAB/GiD figures
  poo.poisson2d_plot_flux.qx='yes';
  poo.poisson2d_plot_flux.qy='yes'; 
  poo.poisson2d_plot_flux.qnorm='yes';   % norm of the flux
  % colorbar limits
  poo.poisson2d_plot_flux.clim.qx=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  poo.poisson2d_plot_flux.clim.qy=[]; % leave an empty [] for matlab's default value, or specify [min max] values.   
  poo.poisson2d_plot_flux.clim.qnorm=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  
  % plotting of gradients to MATLAB/GiD figures
  poo.poisson2d_plot_grad.dx='no';
  poo.poisson2d_plot_grad.dy='no';   
  poo.poisson2d_plot_grad.dnorm='yes';   % norm of the gradient
  % colorbar limits
  poo.poisson2d_plot_grad.clim.dx=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  poo.poisson2d_plot_grad.clim.dy=[]; % leave an empty [] for matlab's default value, or specify [min max] values.  
  poo.poisson2d_plot_grad.clim.dnorm=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  
  %% LINEARELASTICITY MODULE
  
  % options for plotting deformed domain to MATLAB figures
  poo.linelast2d_plot_deformed_domain='yes';
  poo.linelast2d_scale_for_plotting_deformed_domain=1; % a number > 1 will scale the deformed domain when plotting to MATLAB figures

  % plotting of main variables to MATLAB/GiD/VTK figures
  poo.linelast2d_plot_displacement.ux='yes';  
  poo.linelast2d_plot_displacement.uy='yes';  
  poo.linelast2d_plot_displacement.unorm='yes';  % norm of the displacement
  % colorbar limits
  poo.linelast2d_plot_displacement.clim.ux=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  poo.linelast2d_plot_displacement.clim.uy=[]; % leave an empty [] for matlab's default value, or specify [min max] values.   
  poo.linelast2d_plot_displacement.clim.unorm=[]; % leave an empty [] for matlab's default value, or specify [min max] values.   
  
  % plotting of stresses to MATLAB/GiD figures
  poo.linelast2d_plot_stress.s11='yes';
  poo.linelast2d_plot_stress.s12='no';
  poo.linelast2d_plot_stress.s22='no';
  poo.linelast2d_plot_stress.s33='no';
  poo.linelast2d_plot_stress.s1='no';
  poo.linelast2d_plot_stress.s2='no'; 
  poo.linelast2d_plot_stress.s3='no';
  poo.linelast2d_plot_stress.vm='yes';
  poo.linelast2d_plot_stress.p='yes';
  % colorbar limits
  poo.linelast2d_plot_stress.clim.s11=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  poo.linelast2d_plot_stress.clim.s12=[]; % leave an empty [] for matlab's default value, or specify [min max] values.  
  poo.linelast2d_plot_stress.clim.s22=[]; % leave an empty [] for matlab's default value, or specify [min max] values.     
  poo.linelast2d_plot_stress.clim.s33=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  poo.linelast2d_plot_stress.clim.s1=[]; % leave an empty [] for matlab's default value, or specify [min max] values.  
  poo.linelast2d_plot_stress.clim.s2=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  poo.linelast2d_plot_stress.clim.s3=[]; % leave an empty [] for matlab's default value, or specify [min max] values.   
  poo.linelast2d_plot_stress.clim.vm=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  poo.linelast2d_plot_stress.clim.p=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  
  % plotting of strains to MATLAB/GiD figures
  poo.linelast2d_plot_strain.e11='yes';
  poo.linelast2d_plot_strain.e12='no';
  poo.linelast2d_plot_strain.e22='no';
  poo.linelast2d_plot_strain.e33='no';
  poo.linelast2d_plot_strain.e1='no';
  poo.linelast2d_plot_strain.e2='no';  
  poo.linelast2d_plot_strain.e3='no'; 
  % colorbar limits
  poo.linelast2d_plot_strain.clim.e11=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  poo.linelast2d_plot_strain.clim.e12=[]; % leave an empty [] for matlab's default value, or specify [min max] values.  
  poo.linelast2d_plot_strain.clim.e22=[]; % leave an empty [] for matlab's default value, or specify [min max] values.    
  poo.linelast2d_plot_strain.clim.e33=[]; % leave an empty [] for matlab's default value, or specify [min max] values.  
  poo.linelast2d_plot_strain.clim.e1=[]; % leave an empty [] for matlab's default value, or specify [min max] values.   
  poo.linelast2d_plot_strain.clim.e2=[]; % leave an empty [] for matlab's default value, or specify [min max] values.   
  poo.linelast2d_plot_strain.clim.e3=[]; % leave an empty [] for matlab's default value, or specify [min max] values.
  
end

