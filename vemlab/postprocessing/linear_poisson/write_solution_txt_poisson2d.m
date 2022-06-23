function write_solution_txt_poisson2d(domainMesh,scalar_sol,flux,grad,config)

  if strcmp(config.vemlab_method,'VEM2D')||strcmp(config.vemlab_method,'FEM2DT3')
    write_solution_txt_VEM2D_FEM2DT3_poisson2d(domainMesh,scalar_sol,flux,...
                                                grad,config);
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    write_solution_txt_FEM2DQ4_poisson2d(domainMesh,scalar_sol,flux,...
                                          grad,config);    
  else
    throw_error('Error in write_solution_txt_poisson2d.m: vemlab_method')
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:          write_solution_txt_VEM2D_FEM2DT3_poisson2d
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Write "method" solutions to a text file, where "method" can be VEM2D or 
% FEM2DT3. Both methods give constant flux and gradient in the element.
%
% Usage
% =====
% write_solution_txt_VEM2D_FEM2DT3_poisson2d(domainMesh,scalar_sol,flux,...
%                                            grad,config)
%
% Input
% =====
% domainMesh    : structure containing the polygonal mesh information
% scalar_sol    : nodal scalar field solution
% flux          : structure storing flux components
% grad          : structure storing gradient components
% config        : structure storing VEMLab configuration options and behavior
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
% Mar. 17, 2018: first realease (by A. Ortiz-Bernardin)
% Feb. 1, 2020: add an input array variable called triangles_per_polygon, which
%               is used to fix an error in the plotting of VEM fluxes and gradients 
%               into a text file stage (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_solution_txt_VEM2D_FEM2DT3_poisson2d(domainMesh,scalar_sol,...
                                                    flux,grad,config)
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
  fprintf(fid,'*                   Poisson Problem                    *\n');    
  fprintf(fid,'********************************************************\n\n');   
  fprintf(fid,'Mesh filename: %s\n',config.mesh_filename);   
  fprintf(fid,'Method: %s\n',config.vemlab_method);   
  fprintf(fid,'Date: %d-%d-%d\n',time_date(3),time_date(2),time_date(1));
  fprintf(fid,'Time: %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5));  
  fprintf(fid,'Nodal Scalar Field: \n');
  fprintf(fid,' Node      Coords         u \n');
  for i = 1:num_nodes
    fprintf(fid,'%d %.8f %.8f %.8f\n', ...
            i,domainMesh.coords(i,1),domainMesh.coords(i,2),scalar_sol(i));
  end
  % write fluxes and gradients only if flux is not empty  
  if ~isempty(flux) 
    fluxNorm=sqrt((flux.qx).*(flux.qx)+(flux.qy).*(flux.qy));
    gradNorm=sqrt((grad.dx).*(grad.dx)+(grad.dy).*(grad.dy));
    fprintf(fid,'Gradient and Flux: \n'); 
    fprintf(fid,'(quantities are constant in the element)\n\n');   
    numel=length(domainMesh.connect);
    for elem_i=1:numel
      fprintf(fid,'Element %d\n',elem_i);  
      fprintf(fid,'Grad_x      Grad_y      ||Grad||      Flux_x       Flux_y      ||Flux|| \n');
      fprintf(fid,'%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n', ...
              grad.dx(elem_i),grad.dy(elem_i),gradNorm(elem_i),...
              flux.qx(elem_i),flux.qy(elem_i),fluxNorm(elem_i));          
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
% FUNCTION:                 write_solution_txt_FEM2DQ4_poisson2d 
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
% write_solution_txt_FEM2DQ4_poisson2d(domainMesh,scalar_sol,flux,...
%                                      grad,config)
%
% Input
% =====
% domainMesh    : structure containing the polygonal mesh information
% scalar_sol    : nodal scalar field solution
% flux          : structure storing flux components
% grad          : structure storing gradient components
% config        : structure storing VEMLab configuration options and behavior
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
% Mar. 17, 2018: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_solution_txt_FEM2DQ4_poisson2d(domainMesh,scalar_sol,flux,grad,config)
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
  fprintf(fid,'*                   Poisson Problem                    *\n');    
  fprintf(fid,'********************************************************\n\n');   
  fprintf(fid,'Mesh filename: %s\n',config.mesh_filename);   
  fprintf(fid,'Method: %s\n',config.vemlab_method);   
  fprintf(fid,'Date: %d-%d-%d\n',time_date(3),time_date(2),time_date(1));
  fprintf(fid,'Time: %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5));  
  fprintf(fid,'Nodal Displacements: \n');
  fprintf(fid,' Node      Coords         u1       u2 \n');
  for i = 1:num_nodes
    fprintf(fid,'%d %.8f %.8f %.8f\n', ...
            i,domainMesh.coords(i,1),domainMesh.coords(i,2),scalar_sol(i));
  end
  % write fluxes and gradients only if flux is not empty
  if ~isempty(flux)   
    fluxNorm=sqrt((flux.qx).*(flux.qx)+(flux.qy).*(flux.qy));
    gradNorm=sqrt((grad.dx).*(grad.dx)+(grad.dy).*(grad.dy));
    fprintf(fid,'Gradient and Flux: \n'); 
    fprintf(fid,'(at Gauss points --- 2x2 quadrature rule)\n\n');   
    numel=length(domainMesh.connect);
    num_elem_gp=4; % 2x2 rule
    gp=0;
    for elem_i=1:numel
      elem_gp=0;    
      fprintf(fid,'Element %d\n',elem_i);  
      fprintf(fid,'GP#   Grad_x      Grad_y      ||Grad||      Flux_x       Flux_y      ||Flux|| \n');
      for gpoints=1:num_elem_gp
        gp=gp+1;
        elem_gp=elem_gp+1;      
        fprintf(fid,'%d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n',...
                elem_gp,grad.dx(gp),grad.dy(gp),gradNorm(gp),...
                flux.qx(gp),flux.qy(gp),fluxNorm(gp));   
      end
    end
  end
  fprintf('Check txt output files in folder: %s\n',...
           config.txt_output_folder_location);   
  fclose(fid); 
end

