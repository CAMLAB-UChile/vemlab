function write_solution_VTK_linelast2d(domainMesh,displacements,config)
  % plot to a VTK file   
  writeVTK(domainMesh.coords,domainMesh.connect,displacements,config);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:                            writeVTK 
%
% Created by : R. Silva-Valenzuela, rsilvav@outlook.com
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Write solutions to a VTK [1] output file.
%
% Usage
% =====
% writeVTK(coord,conec,uxy,config)
%
% Input
% =====
% coord       : nodal coordinates
% conec       : nodal connectivity
% uxy         : nodal displacement solution
% config      : structure storing VEMLab configuration options and behavior
%
% Output
% ======
%
%-------------------------------------------------------------------------------
% References 
% ==========
% [1] The Visualization Toolkit (VTK), https://www.vtk.org/
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Apr. 9, 2018: first realease (by R. Silva-Valenzuela)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writeVTK(coord,conec,uxy,config)

  fprintf('\n'); 
  fprintf('Writing %s solution to a VTK file...\n',config.vemlab_method); 
  output_filename_results=strcat(config.mesh_filename,'.vtk');
  outfile_results=[config.VTK_output_folder_location,output_filename_results];
  fid=fopen(outfile_results,'wt');   
  fprintf(fid,'# vtk DataFile Version 1.0 \n');
  fprintf(fid,'PolyMesh \n');
  fprintf(fid,'ASCII \n \n');
  fprintf(fid,'DATASET POLYDATA \n \n');
  fprintf(fid,'POINTS ');
  fprintf(fid,'%d', length(coord));
  fprintf(fid,' float \n');
 
  for i=1:length(coord)
    fprintf(fid,'%e ',coord(i,1),coord(i,2), 0.0); 
    fprintf(fid,' \n ');
  end
  fprintf(fid,' \n ');

  c=0;
  for i=1:length(conec)
     c = c + length(conec{i,1}); %conec(i,1)
  end
  c = c + length(conec);

  fprintf(fid,'POLYGONS ');
  fprintf(fid,'%d', length(conec));
  fprintf(fid,' ');
  fprintf(fid,'%d', c);
  fprintf(fid,' \n ');

  for i=1:length(conec)
    fprintf(fid,'%d', length(conec{i,1}));
    fprintf(fid,' ');
    for j=1:length(conec{i,1})
    fprintf(fid,'%d ', conec{i,1}(j)-1);
    fprintf(fid,' ');
    end     
    fprintf(fid,' \n ');
  end
  fprintf(fid,' \n ');

  fprintf(fid,'POINT_DATA ');
  fprintf(fid,'%d', length(coord));
  fprintf(fid,' \n ');

  fprintf(fid,'VECTORS Displacement float \n ');
  for i=1:length(coord)
    fprintf(fid,'%e ',uxy(2*i-1),uxy(2*i), 0.0); 
    fprintf(fid,' \n ');
  end
  fprintf(fid,' \n ');
  
  fprintf('Check VTK output files in folder: %s\n',...
           config.VTK_output_folder_location);   
  fclose(fid);  
  
end

