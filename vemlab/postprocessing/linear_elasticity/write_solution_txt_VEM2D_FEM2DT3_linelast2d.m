%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:          write_solution_txt_VEM2D_FEM2DT3_linelast2d
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Write "method" solutions to a text file, where "method" can be VEM2D or 
% FEM2DT3. Both methods give constant stresses and strains in the element.
%
% Usage
% =====
% write_solution_txt_VEM2D_FEM2DT3_linelast2d(mesh,displacements,stresses,...
%                                             strains,config)
%
% Input
% =====
% mesh          : structure containing mesh data (coords,connect,etc.)
% displacements : nodal displacement solution
% stresses    : structure storing stress tensor components and vMises stress
% strains     : structure storing strain tensor components
% config      : structure storing VEMLab configuration options and behavior
%
% Output
% ======
%
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

function write_solution_txt_VEM2D_FEM2DT3_linelast2d(mesh,displacements,stresses,strains,config)
  fprintf('\n'); 
  fprintf('Writing %s solution to a text file...\n',config.vemlab_method); 
  % output file
  output_filename=strcat('out_',config.mesh_filename);
  outfile=[config.txt_output_folder_location,output_filename];
  fid=fopen(outfile,'w');  
  time_date=clock;  
  % write solution
  num_nodes=size(mesh.coords,1);   
  fprintf(fid,'********************************************************\n');  
  fprintf(fid,'*                      VEMLab %s                      *\n',config.vemlab_version); 
  fprintf(fid,'*      http://camlab.cl/research/software/vemlab       *\n');  
  fprintf(fid,'*------------------------------------------------------*\n');    
  fprintf(fid,'*               Linear Elastostatics                   *\n');    
  fprintf(fid,'********************************************************\n\n');   
  fprintf(fid,'Mesh filename: %s\n',config.mesh_filename);   
  fprintf(fid,'Method: %s\n',config.vemlab_method);   
  fprintf(fid,'Date: %d-%d-%d\n',time_date(3),time_date(2),time_date(1));
  fprintf(fid,'Time: %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5));  
  fprintf(fid,'Nodal Displacements: \n');
  fprintf(fid,' Node      Coords         u1       u2 \n');
  for i = 1:num_nodes
    fprintf(fid,'%d %.8f %.8f %.8f %.8f\n', ...
            i,mesh.coords(i,1),mesh.coords(i,2),displacements(2*i-1),displacements(2*i));
  end
  fprintf(fid,'Strains, Stresses, Principal Strains and Principal Stresses: \n'); 
  fprintf(fid,'(quantities are constant in the element)\n\n');   
  numel=length(mesh.connect);
  for elem_i=1:numel
    fprintf(fid,'Element %d\n',elem_i);  
    fprintf(fid,'e_11      e_22      e_33     e_12      s_11       s_22       s_33      s_12      VM      e_1      e_2      e_3      s_1      s_2      s_3 \n');
    fprintf(fid,'%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n', ...
            strains.e11(elem_i),strains.e22(elem_i),strains.e33(elem_i),...
            strains.e12(elem_i),stresses.s11(elem_i),stresses.s22(elem_i),...
            stresses.s33(elem_i),stresses.s12(elem_i),stresses.vm(elem_i),...
            strains.e1(elem_i),strains.e2(elem_i),strains.e3(elem_i),...
            stresses.s1(elem_i),stresses.s2(elem_i),stresses.s3(elem_i));   
  end
  fprintf('Check txt output files in folder: %s\n',...
           config.txt_output_folder_location);   
  fclose(fid); 
end