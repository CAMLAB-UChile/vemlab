%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                     config_vemlab
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Configure VEMLab program options and behavior
%
% Usage
% =====
% config = config_vemlab(opsystem,vemlab_root_dir,mesh_filename,problem,...
%                        plot_mesh_over_results,write_solutions_to_text_file,...
%                        write_solutions_to_GiD_file,vemlab_module,...
%                        vemlab_method)
%
% Input
% =====
% opsystem               : machine's operating system
% vemlab_root_dir        : root directory for VEMLab
% mesh_file_name         : file that contains the mesh information
% plot_mesh_over_results : 'yes' or 'no'
% write_solutions_to_text_file : 'yes' or 'no'
% write_solutions_to_GiD_file : 'yes' or 'no' 
% vemlab_module          : 'LinearElastostatics' or 'Poisson'
% vemlab_method          : 'VEM2D', 'FEM2DT3' or 'FEM2DQ4'
%
% Output
% ======
% config : structure storing the input variables
%-------------------------------------------------------------------------------
% References 
% ==========
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Oct. 19, 2018: add variable "config.create_matlab_contour_plots" to control
%                whether matlab contour plots must be created or not 
%                (by A. Ortiz-Bernardin)
% May 13, 2018: add variable "config.number_of_gauss_points_per_axis_FEM2DQ4"
%               to set the number of Gauss points to integrate the FEM2DQ4 
%               stiffness matrix and body force vector (by A. Ortiz-Bernardin)
% Mar. 22, 2018: add options for plotting purposes (by A. Ortiz-Bernardin)
% Mar. 17, 2018: add Poisson module (by A. Ortiz-Bernardin)
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function config = config_vemlab(opsystem,vemlab_root_dir,mesh_filename,...
                                vemlab_module,vemlab_method)
  % VEMLab version
  config.vemlab_version='2.2.2';
  
  % program options
  [create_matlab_contour_plots,plot_mesh,plot_mesh_over_results,write_solutions_to_text_file,...
   write_solutions_to_GiD_file,write_solutions_to_VTK_file,poisson2d_plot_scalar_field,...
   poisson2d_plot_flux,poisson2d_plot_grad,linelast2d_plot_displacement,...
   linelast2d_plot_stress,linelast2d_plot_strain,linelast2d_plot_deformed_domain,...
   linelast2d_scale_for_plotting_deformed_domain]=plot_and_output_options; 
  
  config.opsystem=opsystem;
  config.vemlab_root_dir=vemlab_root_dir;  
  config.mesh_filename=mesh_filename;
  config.create_matlab_contour_plots=create_matlab_contour_plots;   
  config.plot_mesh=plot_mesh;  
  config.plot_mesh_over_results=plot_mesh_over_results;   
  config.write_solutions_to_text_file=write_solutions_to_text_file;  
  config.write_solutions_to_GiD_file=write_solutions_to_GiD_file; 
  config.write_solutions_to_VTK_file=write_solutions_to_VTK_file;  
  config.poisson2d_plot_scalar_field=poisson2d_plot_scalar_field;  
  config.poisson2d_plot_flux=poisson2d_plot_flux;
  config.poisson2d_plot_grad=poisson2d_plot_grad;
  config.linelast2d_plot_displacement=linelast2d_plot_displacement;  
  config.linelast2d_plot_stress=linelast2d_plot_stress;
  config.linelast2d_plot_strain=linelast2d_plot_strain;
  config.vemlab_module=vemlab_module;  
  config.vemlab_method=vemlab_method;
  
  config.linelast2d_plot_deformed_domain=linelast2d_plot_deformed_domain;
  config.linelast2d_scale_for_plotting_deformed_domain=linelast2d_scale_for_plotting_deformed_domain;
  
  config.number_of_gauss_points_per_axis_FEM2DQ4=2; % for Poisson and 
                                                    % linear elastostatic problems, 
                                                    % this should be set to 2 as
                                                    % a minimum. 1 is for
                                                    % underintegration
    
  logical_cond_poisson2d=strcmp(poisson2d_plot_flux.qx,'yes')||...
                         strcmp(poisson2d_plot_flux.qy,'yes') ||...
                         strcmp(poisson2d_plot_flux.qnorm,'yes') ||...
                         strcmp(poisson2d_plot_grad.dx,'yes') ||...
                         strcmp(poisson2d_plot_grad.dy,'yes') ||...
                         strcmp(poisson2d_plot_grad.dnorm,'yes');
  if logical_cond_poisson2d
    config.poisson2d_plot_flux_and_gradient='yes';
  else
    config.poisson2d_plot_flux_and_gradient='no';
  end
  
  logical_cond_linelast2d=strcmp(linelast2d_plot_stress.s11,'yes') ||...
                          strcmp(linelast2d_plot_stress.s12,'yes') ||...
                          strcmp(linelast2d_plot_stress.s22,'yes') ||...                          
                          strcmp(linelast2d_plot_stress.s33,'yes') ||...
                          strcmp(linelast2d_plot_stress.s1,'yes') ||...
                          strcmp(linelast2d_plot_stress.s2,'yes') ||...
                          strcmp(linelast2d_plot_stress.s3,'yes') ||...      
                          strcmp(linelast2d_plot_stress.vm,'yes') ||... 
                          strcmp(linelast2d_plot_strain.e11,'yes') ||...   
                          strcmp(linelast2d_plot_strain.e12,'yes') ||...
                          strcmp(linelast2d_plot_strain.e22,'yes') ||...     
                          strcmp(linelast2d_plot_strain.e33,'yes') ||...  
                          strcmp(linelast2d_plot_strain.e1,'yes') ||...  
                          strcmp(linelast2d_plot_strain.e2,'yes') ||... 
                          strcmp(linelast2d_plot_strain.e3,'yes');                         
  if logical_cond_linelast2d
    config.linelast2d_plot_stress_and_strain='yes';
  else
    config.linelast2d_plot_stress_and_strain='no';
  end  
  
  % input/output folders
  is_Windows = strcmp(opsystem,'PCWIN') || strcmp(opsystem,'PCWIN64');  
  is_Linux = strcmp(opsystem,'GLNX86') || strcmp(opsystem,'GLNXA64');  
  if is_Windows
    config.txt_output_folder_location=[vemlab_root_dir,'\test\output_files\txt\'];
    config.GiD_output_folder_location=[vemlab_root_dir,'\test\output_files\GiD\'];
    config.VTK_output_folder_location=[vemlab_root_dir,'\test\output_files\VTK\'];      
    config.mesh_folder_location=[vemlab_root_dir,'\test\mesh_files\'];
  elseif is_Linux
    config.txt_output_folder_location=[vemlab_root_dir,'/test/output_files/txt/'];
    config.GiD_output_folder_location=[vemlab_root_dir,'/test/output_files/GiD/'];
    config.VTK_output_folder_location=[vemlab_root_dir,'/test/output_files/VTK/'];     
    config.mesh_folder_location=[vemlab_root_dir,'/test/mesh_files/'];  
  end 
  
  % some global checks
  logical_cond=strcmp(config.vemlab_module,'LinearElastostatics')||strcmp(config.vemlab_module,'Poisson');  
  if ~logical_cond
    throw_error('Error in config_vemlab.m: vemlab_module\n');
  end
  
  logical_cond=strcmp(config.vemlab_method,'VEM2D')||strcmp(config.vemlab_method,'FEM2DT3');
  logical_cond=logical_cond||strcmp(config.vemlab_method,'FEM2DQ4');
  if ~logical_cond
    throw_error('Error in config_vemlab.m: vemlab_method\n');
  end  
  
end

