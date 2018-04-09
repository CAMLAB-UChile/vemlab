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
%                        plot_mesh_over_results,vemlab_module,vemlab_method)
%
% Input
% =====
% opsystem               : machine's operating system
% vemlab_root_dir        : root directory for VEMLab
% mesh_file_name         : file that contains the mesh information
% plot_mesh_over_results : 'yes' or 'no'
% vemlab_module          : 'LinearElastostatics'
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
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function config = config_vemlab(opsystem,vemlab_root_dir,mesh_filename,...
                                plot_mesh_over_results,vemlab_module,...
                                vemlab_method)
  % VEMLab version
  config.vemlab_version='1.0';
  
  % program options
  config.opsystem=opsystem;
  config.vemlab_root_dir=vemlab_root_dir;  
  config.mesh_filename=mesh_filename; 
  config.plot_mesh_over_results=plot_mesh_over_results;       
  config.vemlab_module=vemlab_module;  
  config.vemlab_method=vemlab_method;
  
  % input/output folders
  is_Windows = strcmp(opsystem,'PCWIN') || strcmp(opsystem,'PCWIN64');  
  is_Linux = strcmp(opsystem,'GLNX86') || strcmp(opsystem,'GLNXA64');  
  if is_Windows
    config.txt_output_folder_location=[vemlab_root_dir,'\test\output_files\txt\'];
    config.GiD_output_folder_location=[vemlab_root_dir,'\test\output_files\GiD\'];    
    config.mesh_folder_location=[vemlab_root_dir,'\test\mesh_files\'];
  elseif is_Linux
    config.txt_output_folder_location=[vemlab_root_dir,'/test/output_files/txt/'];
    config.GiD_output_folder_location=[vemlab_root_dir,'/test/output_files/GiD/']; 
    config.mesh_folder_location=[vemlab_root_dir,'/test/mesh_files/'];  
  end 
  
  % some global checks
  if ~strcmp(config.vemlab_module,'LinearElastostatics')
    throw_error('Error in config_vemlab.m: vemlab_module\n');
  end
  
  logical_cond=strcmp(config.vemlab_method,'VEM2D')||strcmp(config.vemlab_method,'FEM2DT3');
  logical_cond=logical_cond||strcmp(config.vemlab_method,'FEM2DQ4');
  if ~logical_cond
    throw_error('Error in config_vemlab.m: vemlab_method\n');
  end  
  
end

