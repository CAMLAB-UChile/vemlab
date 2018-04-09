function [plot_mesh,plot_mesh_over_results,write_solutions_to_text_file,...
          write_solutions_to_GiD_file,write_solutions_to_VTK_file,poisson2d_plot_scalar_field,...
          poisson2d_plot_flux,poisson2d_plot_grad,linelast2d_plot_displacement,...
          linelast2d_plot_stress,linelast2d_plot_strain] = plot_and_output_options
  
  % INSTRUCTIONS: SET VARIABLES WITH 'yes' or 'no' ONLY.
  
  %% GENERAL
  
  plot_mesh='yes';
  plot_mesh_over_results='no';
  write_solutions_to_text_file='yes';
  write_solutions_to_GiD_file='yes';
  write_solutions_to_VTK_file='yes';  
  
  %% POISSON MODULE
  
  % plotting of main variables
  poisson2d_plot_scalar_field.u='yes';
  
  % plotting of fluxes 
  poisson2d_plot_flux.qx='no';
  poisson2d_plot_flux.qy='no'; 
  poisson2d_plot_flux.qnorm='yes';   % norm of the flux
  
  % plotting of gradients
  poisson2d_plot_grad.dx='no';
  poisson2d_plot_grad.dy='no';   
  poisson2d_plot_grad.dnorm='yes';   % norm of the gradient
  
  %% LINELAST MODULE

  % plotting of main variables
  linelast2d_plot_displacement.ux='yes';  
  linelast2d_plot_displacement.uy='yes';  
  linelast2d_plot_displacement.unorm='yes';  % norm of the displacement
  
  % plotting of stresses
  linelast2d_plot_stress.s11='yes';
  linelast2d_plot_stress.s12='no';
  linelast2d_plot_stress.s22='no';
  linelast2d_plot_stress.s33='no';
  linelast2d_plot_stress.s1='no';
  linelast2d_plot_stress.s2='no';  
  linelast2d_plot_stress.s3='no';  
  linelast2d_plot_stress.vm='no';  
  
  % plotting of strains
  linelast2d_plot_strain.e11='no';
  linelast2d_plot_strain.e12='no';
  linelast2d_plot_strain.e22='no';
  linelast2d_plot_strain.e33='no';
  linelast2d_plot_strain.e1='no';
  linelast2d_plot_strain.e2='no';  
  linelast2d_plot_strain.e3='no';   
  
end

