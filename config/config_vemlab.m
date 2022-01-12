%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                      VEMLab
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
%                        vemlab_method,vemlab_solver)
%
% Input
% =====
% opsystem               : machine's operating system
% vemlab_root_dir        : root directory for VEMLab
% mesh_file_name         : file that contains the mesh information
% plot_mesh_over_results : 'yes' or 'no'
% write_solutions_to_text_file : 'yes' or 'no'
% write_solutions_to_GiD_file : 'yes' or 'no' 
% vemlab_module          : 'LinearElastostatics', 'LinearElasticFractureXVEM' or 'Poisson'
% vemlab_method          : 'VEM2D', 'FEM2DT3' or 'FEM2DQ4'
% vemlab_solver          : 'sparse' or 'dense'
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
% Jan. 04, 2021: plotting functionalities improved
% May 7, 2020: - add sparse solver for FEM in LinearElastostatics module
%              - add sparse solver for VEM and FEM in Poisson module
%              - add NonlinearPoisson module
% Feb. 10, 2020: add sparse solver for VEM in LinearElastostatics module
% Feb. 2, 2020: add some control variables for mesh plotting; add option to
%               write a result file to be read a Convex Polygon Packing (CPP) program
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
                                vemlab_module,vemlab_method,vemlab_solver,...
                                stability_type)
  % --------------------------------------------------------------------------------------
  % VEMLab version
  % --------------------------------------------------------------------------------------
  config.vemlab_version='2.4';
  
  % ---------------------------
  % plot and output options
  % ---------------------------
  poo=plot_and_output_options; 
  
  % --------------------------------------------------------------------------------------
  % Common options
  % --------------------------------------------------------------------------------------
  config.opsystem=opsystem;
  config.vemlab_root_dir=vemlab_root_dir;  
  config.mesh_filename=mesh_filename;
  config.create_matlab_contour_plots=poo.create_matlab_contour_plots;  
  config.create_matlab_exact_contour_plots=poo.create_matlab_exact_contour_plots;
  config.front_matlab_plot_view_orientation=poo.front_matlab_plot_view_orientation;
  config.plot_mesh_over_results=poo.plot_mesh_over_results;
  config.plot_figure_axis=poo.plot_figure_axis;
  config.matlab_colormap=poo.matlab_colormap;
  config.colormap_number_of_colors=poo.colormap_number_of_colors;
  config.colorbar_tick_label_notation=poo.colorbar_tick_label_notation;
  config.colorbar_TickLength=poo.colorbar_TickLength;
  config.colorbar_LineWidth=poo.colorbar_LineWidth;
  config.axis_xtick=poo.axis_xtick;
  config.axis_ytick=poo.axis_ytick;
  config.axis_ztick=poo.axis_ztick;  
  config.figure_XMinorTick=poo.figure_XMinorTick;
  config.figure_YMinorTick=poo.figure_YMinorTick;
  config.figure_ZMinorTick=poo.figure_ZMinorTick;   
  config.figure_XGrid=poo.figure_XGrid;
  config.figure_YGrid=poo.figure_YGrid;
  config.figure_ZGrid=poo.figure_ZGrid;
  config.figure_XMinorGrid=poo.figure_XMinorGrid;
  config.figure_YMinorGrid=poo.figure_YMinorGrid;    
  config.figure_ZMinorGrid=poo.figure_ZMinorGrid;
  config.figure_GridColor=poo.figure_GridColor;
  config.figure_MinorGridColor=poo.figure_MinorGridColor;
  config.print_figures=poo.print_figures;
  config.print_exact_figures=poo.print_exact_figures;
  config.print_figures_resolution=poo.print_figures_resolution; 
  config.save_matlab_figures=poo.save_matlab_figures;
  config.save_exact_matlab_figures=poo.save_exact_matlab_figures;
  config.plot_vemlab_logo=poo.plot_vemlab_logo;
  config.plot_mesh=poo.plot_mesh; 
  config.plot_mesh_linewidth=poo.plot_mesh_linewidth;  
  config.plot_mesh_nodes=poo.plot_mesh_nodes;
  config.plot_mesh_nodesize=poo.plot_mesh_nodesize;
  config.plot_mesh_axis=poo.plot_mesh_axis;     
  config.write_solutions_to_text_file=poo.write_solutions_to_text_file;  
  config.write_solutions_to_GiD_file=poo.write_solutions_to_GiD_file; 
  config.write_solutions_to_VTK_file=poo.write_solutions_to_VTK_file;
  config.write_solutions_to_CPP_file=poo.write_solutions_to_CPP_file;
  config.vemlab_module=vemlab_module;  
  config.vemlab_method=vemlab_method;
  config.vemlab_solver=vemlab_solver;
  config.stability_type=stability_type;
  config.number_of_gauss_points_per_axis_FEM2DQ4=2; % for Poisson and 
                                                    % linear elastostatic problems, 
                                                    % this should be set to 2 as
                                                    % a minimum. 1 is for
                                                    % underintegration
    
  % --------------------------------------------------------------------------------------                                                
  % Poisson module options
  % --------------------------------------------------------------------------------------
  config.poisson2d_plot_scalar_field=poo.poisson2d_plot_scalar_field;  
  config.poisson2d_plot_flux=poo.poisson2d_plot_flux;
  config.poisson2d_plot_grad=poo.poisson2d_plot_grad;
  
  logical_cond_poisson2d=strcmp(poo.poisson2d_plot_flux.qx,'yes')||...
                         strcmp(poo.poisson2d_plot_flux.qy,'yes') ||...
                         strcmp(poo.poisson2d_plot_flux.qnorm,'yes') ||...
                         strcmp(poo.poisson2d_plot_grad.dx,'yes') ||...
                         strcmp(poo.poisson2d_plot_grad.dy,'yes') ||...
                         strcmp(poo.poisson2d_plot_grad.dnorm,'yes');
  if logical_cond_poisson2d
    config.poisson2d_plot_flux_and_gradient='yes';
  else
    config.poisson2d_plot_flux_and_gradient='no';
  end
  
  % --------------------------------------------------------------------------------------
  % LinearElasticity module options
  % --------------------------------------------------------------------------------------
  config.linelast2d_plot_displacement=poo.linelast2d_plot_displacement;  
  config.linelast2d_plot_stress=poo.linelast2d_plot_stress;
  config.linelast2d_plot_strain=poo.linelast2d_plot_strain;
  config.linelast2d_plot_deformed_domain=poo.linelast2d_plot_deformed_domain;
  config.linelast2d_scale_for_plotting_deformed_domain=poo.linelast2d_scale_for_plotting_deformed_domain;  
  
  logical_cond_linelast2d=strcmp(poo.linelast2d_plot_stress.s11,'yes') ||...
                          strcmp(poo.linelast2d_plot_stress.s12,'yes') ||...
                          strcmp(poo.linelast2d_plot_stress.s22,'yes') ||...                          
                          strcmp(poo.linelast2d_plot_stress.s33,'yes') ||...
                          strcmp(poo.linelast2d_plot_stress.s1,'yes') ||...
                          strcmp(poo.linelast2d_plot_stress.s2,'yes') ||...
                          strcmp(poo.linelast2d_plot_stress.s3,'yes') ||...      
                          strcmp(poo.linelast2d_plot_stress.vm,'yes') ||... 
                          strcmp(poo.linelast2d_plot_stress.p,'yes') ||...                           
                          strcmp(poo.linelast2d_plot_strain.e11,'yes') ||...   
                          strcmp(poo.linelast2d_plot_strain.e12,'yes') ||...
                          strcmp(poo.linelast2d_plot_strain.e22,'yes') ||...     
                          strcmp(poo.linelast2d_plot_strain.e33,'yes') ||...  
                          strcmp(poo.linelast2d_plot_strain.e1,'yes') ||...  
                          strcmp(poo.linelast2d_plot_strain.e2,'yes') ||... 
                          strcmp(poo.linelast2d_plot_strain.e3,'yes');                         
  if logical_cond_linelast2d
    config.linelast2d_plot_stress_and_strain='yes';
  else
    config.linelast2d_plot_stress_and_strain='no';
  end
 
  % --------------------------------------------------------------------------------------
  % Input/Output folders
  % --------------------------------------------------------------------------------------
  is_Windows = strcmp(opsystem,'PCWIN') || strcmp(opsystem,'PCWIN64');  
  is_Linux = strcmp(opsystem,'GLNX86') || strcmp(opsystem,'GLNXA64');  
  if is_Windows
    config.txt_output_folder_location=[vemlab_root_dir,'\test\output_files\txt\'];
    config.GiD_output_folder_location=[vemlab_root_dir,'\test\output_files\GiD\'];
    config.VTK_output_folder_location=[vemlab_root_dir,'\test\output_files\VTK\'];    
    config.CPP_output_folder_location=[vemlab_root_dir,'\test\output_files\CPP\'];    
    config.matlab_figures_output_folder_location=[vemlab_root_dir,'\test\output_files\matlab_figures\'];     
    config.mesh_folder_location=[vemlab_root_dir,'\test\mesh_files\'];
    config.GiD_T3_mesh_folder_location=[vemlab_root_dir,'\test\mesh_files\GiD_T3_mesh_files\']; 
    config.internal_files_folder_location=[vemlab_root_dir,'\vemlab\internal_files\'];  
  elseif is_Linux
    config.txt_output_folder_location=[vemlab_root_dir,'/test/output_files/txt/'];
    config.GiD_output_folder_location=[vemlab_root_dir,'/test/output_files/GiD/'];
    config.VTK_output_folder_location=[vemlab_root_dir,'/test/output_files/VTK/'];  
    config.CPP_output_folder_location=[vemlab_root_dir,'/test/output_files/CPP/'];  
    config.matlab_figures_output_folder_location=[vemlab_root_dir,'/test/output_files/matlab_figures/'];    
    config.mesh_folder_location=[vemlab_root_dir,'/test/mesh_files/'];  
    config.GiD_T3_mesh_folder_location=[vemlab_root_dir,'/test/mesh_files/GiD_T3_mesh_files/'];
    config.internal_files_folder_location=[vemlab_root_dir,'/vemlab/internal_files/'];
  end 
  
  % --------------------------------------------------------------------------------------
  % Global checks
  % --------------------------------------------------------------------------------------
  logical_cond=strcmp(config.vemlab_module,'LinearElastostatics')||...
               strcmp(config.vemlab_module,'Poisson');
  if ~logical_cond
    throw_error('Error in config_vemlab.m: vemlab_module\n');
  end
  
  logical_cond=strcmp(config.vemlab_method,'VEM2D')||strcmp(config.vemlab_method,'FEM2DT3');
  logical_cond=logical_cond||strcmp(config.vemlab_method,'FEM2DQ4');
  if ~logical_cond
    throw_error('Error in config_vemlab.m: vemlab_method\n');
  end  
  
  % --------------------------------------------------------------------------------------
  % global variables
  % --------------------------------------------------------------------------------------
  config.gid_mesh_file_is_available=false; % for nonlinear problems with steps this ensures 
                                           % that the GiD mesh file is
                                           % written only one time. As soon
                                           % as the file was written this
                                           % variable will be switched to
                                           % 'true' stating the .msh file
                                           % is already available
  config.gid_result_file_is_available=false; % similar to config.gid_mesh_file_is_available
    
  %---------------------------------------------------------------------------------------
  % delete previous output files and figures
  %---------------------------------------------------------------------------------------
  if strcmp(poo.delete_existing_output_figures,'yes')
    cd(config.matlab_figures_output_folder_location);   
    delete *.fig
    delete *.pdf
    cd(config.vemlab_root_dir);
    cd test;
  end
  if strcmp(poo.delete_existing_GiD_output_files,'yes')
    cd(config.GiD_output_folder_location);   
    delete *.msh
    delete *.res
    cd(vemlab_root_dir);
    cd test;    
  end
  if strcmp(poo.delete_existing_txt_output_files,'yes')
    cd(config.txt_output_folder_location);   
    delete *.txt
    cd(vemlab_root_dir);
    cd test;    
  end  
  if strcmp(poo.delete_existing_VTK_output_files,'yes')
    cd(config.VTK_output_folder_location);   
    delete *.vtk
    cd(vemlab_root_dir);
    cd test;    
  end   
  
  %%%---------------------------------------------------------------------------
  %%% WARNING: DO NOT CHANGE THE LINES BELOW. THIS IS EXPERIMENTAL ONLY.
  %%% SEE THE DISCUSION IN FUNCTION "print_figure.m" THAT IS LOCATED IN
  %%% FOLDER "postprocessing".
  config.printer_renderer='opengl'; % 'opengl' or 'painters'
  config.path_to_inkscape='';
  % config.path_to_inkscape='"C:\Program Files\Inkscape\inkscape.exe"';   
  % if not installed leave ''. Notice the double quotation between the single quotation.  
  if isempty(config.path_to_inkscape) && strcmp(config.printer_renderer,'painters')
    throw_warning('Painters renderer needs Inkscape. No path to Inkscape was provided. Switching to opengl renderer.');
    fprintf('\n');
    config.printer_renderer='opengl';
  end  
  %%%----------------------------------------------------------------------------  
    
end

