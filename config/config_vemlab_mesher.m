%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                     config_vemlab_mesher
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Configure VEMLab mesher programs options and behavior
%
% Usage
% =====
% config = config_vemlab(opsystem,vemlab_root_dir,mesh_filename)
%
% Input
% =====
% opsystem               : machine's operating system
% vemlab_root_dir        : root directory for VEMLab
% mesh_file_name         : file that contains the mesh information
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

function config = config_vemlab_mesher(opsystem,vemlab_root_dir,mesh_filename)

  % program options
  config.opsystem=opsystem;  
  config.vemlab_root_dir=vemlab_root_dir;    
  config.mesh_filename=mesh_filename; 
  
  % input/output folders
  is_Windows = strcmp(opsystem,'PCWIN') || strcmp(opsystem,'PCWIN64');
  is_Linux = strcmp(opsystem,'GLNX86') || strcmp(opsystem,'GLNXA64');  
  if is_Windows  
    config.mesh_folder_location=[vemlab_root_dir,'\test\mesh_files\'];
  elseif is_Linux 
    config.mesh_folder_location=[vemlab_root_dir,'/test/mesh_files/'];  
  end 
  
end

