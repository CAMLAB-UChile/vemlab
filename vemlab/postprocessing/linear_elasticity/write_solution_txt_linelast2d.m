function write_solution_txt_linelast2d(domainMesh,displacements,stress,strain,config)

  if strcmp(config.vemlab_method,'VEM2D')||strcmp(config.vemlab_method,'FEM2DT3')
    write_solution_txt_VEM2D_FEM2DT3_linelast2d(domainMesh,displacements,stress,...
                                                strain,config);
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    write_solution_txt_FEM2DQ4_linelast2d(domainMesh,displacements,stress,...
                                          strain,config);    
  else
    throw_error('Error in write_solution_txt_linelast2d.m: vemlab_method')
  end

end

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
% write_solution_txt_VEM2D_FEM2DT3_linelast2d(domainMesh,displacements,stresses,...
%                                             strains,config)
%
% Input
% =====
% domainMesh : structure containing mesh data (coords,connect,etc.)
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
% Feb. 1, 2020: add an input array variable called triangles_per_polygon, which
%               is used to fix an error in the plotting of VEM stresses and strains 
%               into a text file stage (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_solution_txt_VEM2D_FEM2DT3_linelast2d(domainMesh,displacements,...
                                                     stresses,strains,config)
  fprintf('\n'); 
  fprintf('Writing %s solution to a text file...\n',config.vemlab_method); 
  % output file
  output_filename=strcat('out_',config.mesh_filename);
  outfile=[config.txt_output_folder_location,output_filename];
  fid=fopen(outfile,'w');  
  time_date=clock;  
  % write solution
  num_nodes=size(domainMesh.coords,1);   
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
            i,domainMesh.coords(i,1),domainMesh.coords(i,2),displacements(2*i-1),displacements(2*i));
  end
  % write stresses and strains only if stresses is not empty.
  if ~isempty(stresses)
    fprintf(fid,'Strains, Stresses, Principal Strains and Principal Stresses: \n'); 
    fprintf(fid,'(quantities are constant in the element)\n\n');   
    numel=length(domainMesh.connect);
    k = 1;
    for elem_i=1:numel
      fprintf(fid,'Element %d\n',elem_i);  
      fprintf(fid,'e_11      e_22      e_33     e_12      s_11       s_22       s_33      s_12      VM      p      e_1      e_2      e_3      s_1      s_2      s_3 \n');
      fprintf(fid,'%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n', ...
              strains.e11(k),strains.e22(k),strains.e33(k),...
              strains.e12(k),stresses.s11(k),stresses.s22(k),...
              stresses.s33(k),stresses.s12(k),stresses.vm(k),stresses.p(k),...
              strains.e1(k),strains.e2(k),strains.e3(k),...
              stresses.s1(k),stresses.s2(k),stresses.s3(k));   
      k = k + 1;
    end
  end
  fprintf('Check txt output files in folder: %s\n',...
           config.txt_output_folder_location);   
  fclose(fid); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                 write_solution_txt_FEM2DQ4_linelast 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Write "FEM2DQ4" solutions to a text file
%
% Usage
% =====
% write_solution_txt_FEM2DQ4_linelast2d(domainMesh,displacements,stresses,...
%                                       strains,config)
%
% Input
% =====
% domainMesh : structure containing mesh data (coords,connect,etc.)
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

function write_solution_txt_FEM2DQ4_linelast2d(domainMesh,displacements,...
                                               stresses,strains,config)
  fprintf('\n'); 
  fprintf('Writing %s solution to a text file...\n',config.vemlab_method); 
  % output file
  output_filename=strcat('out_',config.mesh_filename);
  outfile=[config.txt_output_folder_location,output_filename];
  fid=fopen(outfile,'w');  
  time_date=clock;  
  % write solution
  num_nodes=size(domainMesh.coords,1);   
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
            i,domainMesh.coords(i,1),domainMesh.coords(i,2),displacements(2*i-1),displacements(2*i));
  end
  % write stresses and strains only if stresses is not empty
  if ~isempty(stresses)  
    fprintf(fid,'Strains, Stresses, Principal Strains and Principal Stresses: \n'); 
    fprintf(fid,'(at Gauss points --- 2x2 quadrature rule)\n\n');   
    numel=length(domainMesh.connect);
    num_elem_gp=4; % 2x2 rule
    gp=0;
    for elem_i=1:numel
      elem_gp=0;
      fprintf(fid,'Element %d\n',elem_i);  
      fprintf(fid,'GP#   e_11      e_22      e_33     e_12      s_11       s_22       s_33      s_12      VM      p      e_1      e_2      e_3      s_1      s_2      s_3 \n');
      for gpoints=1:num_elem_gp
        gp=gp+1;
        elem_gp=elem_gp+1;
        fprintf(fid,'%d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n', ...
                elem_gp,strains.e11(gp),strains.e22(gp),strains.e33(gp),...
                strains.e12(gp),stresses.s11(gp),stresses.s22(gp),...
                stresses.s33(gp),stresses.s12(gp),stresses.vm(gp),stresses.p(gp),...
                strains.e1(gp),strains.e2(gp),strains.e3(gp),...
                stresses.s1(gp),stresses.s2(gp),stresses.s3(gp));   
      end
    end
  end
  fprintf('Check txt output files in folder: %s\n',...
           config.txt_output_folder_location);   
  fclose(fid); 
end
