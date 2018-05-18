function [plot_mesh,plot_mesh_over_results,write_solutions_to_text_file,...
          write_solutions_to_GiD_file,write_solutions_to_VTK_file,poisson2d_plot_scalar_field,...
          poisson2d_plot_flux,poisson2d_plot_grad,linelast2d_plot_displacement,...
          linelast2d_plot_stress,linelast2d_plot_strain] = plot_and_output_options
  
  % INSTRUCTIONS: SET VARIABLES WITH 'yes' or 'no' ONLY.
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
  %       and gradients, it is highly recommended using the GiD output files to
  %       visualize them in "GiD the pre and postprocessor" (www.gidhome.com)
  %
  disp("*** DON'T FORGET TO SETUP THE PLOT AND OUTPUT OPTIONS ***");
  disp(" ");
  disp("To setup them go to file config\plot_and_output_options.m");
  disp(" ");  
  disp("NOTE: Currently, the MATLAB plotting of stresses and strains"); 
  disp("(linear elastostatic problem) or fluxes and gradients (Poisson problem)");
  disp("is very limited and won't recognize holes that might come with");
  disp("the domain geometries. It will work fine with rectangular geometries");
  disp("without holes, but the mesh will need to be very refined to obtain");
  disp("quality colormaps on the whole geometry ... just try a coarse mesh");
  disp("and then a refined mesh to see what is being said.");
  disp(" ");    
  disp("If quality colormap plots are required for stresses, strains, fluxes");
  disp("and gradients, it is higly recommended using the GiD output files to");
  disp("visualize them in GiD the pre and postprocessor (www.gidhome.com)");
  disp(" ");    
  disp("Press any key to continue ...");  
  disp(" ");    
  pause;
  
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

