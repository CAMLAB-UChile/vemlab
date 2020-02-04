function write_solution_CPP_linelast2d(domainMesh,displacements,stress,strain,triangles_per_polygon,config)

  if strcmp(config.vemlab_method,'VEM2D')
    write_solution_CPP_VEM2D_linelast2d(domainMesh,displacements,stress,...
                                        strain,triangles_per_polygon,config);
    write_polymesh_CPP_VEM2D_linelast2d(domainMesh,config);
  else
    throw_warning('*** Warning in write_solution_CPP_linelast2d.m: skipping writing of CPP file. Procedure not available for this vemlab_method ***')
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                  VEMLab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:          write_solution_CPP_VEM2D_FEM2DT3_linelast2d
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Write "method" solutions to a text file to be read in "Convex Polygon Packing (CPP)"
% software, where "method" is VEM2D. This method gives constant stresses and strains 
% in the element.
%
% Usage
% =====
% write_solution_CPP_VEM2D_linelast2d(domainMesh,displacements,stresses,...
%                                     strains,config)
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
% Jan. 25, 2020: first realease (by A. Ortiz-Bernardin)
% Feb. 1, 2020: add an input array variable called triangles_per_polygon, which
%               is used to fix an error in the plotting of VEM stresses and strains 
%               into a text file stage (by A. Ortiz-Bernardin); add function 
%               to write the polygonal mesh to a text file in the format expected
%               by the Convex Polygonal Packing program.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_solution_CPP_VEM2D_linelast2d(domainMesh,displacements,...
                                             stresses,strains,...
                                             triangles_per_polygon,config)
  fprintf('\n'); 
  fprintf('Writing %s solution to a CPP text file...\n',config.vemlab_method); 
  % output file
  output_filename=strcat('outCPP_',config.mesh_filename);
  outfile=[config.CPP_output_folder_location,output_filename];
  fid=fopen(outfile,'w');  
  % write solution
  num_nodes=size(domainMesh.coords,1);   
  fprintf(fid,'ResultName: %s\n','Displacements');   
  fprintf(fid,'NumberOfComponents: %d\n',3);
  fprintf(fid,'ComponentsNames: %s %s %s\n','Ux','Uy','Norm');  
  fprintf(fid,'NumberOfVertices: %d\n',num_nodes);  
  fprintf(fid,' Vertex Ux Uy Norm\n');
  for i = 1:num_nodes
    fprintf(fid,'%d %.16f %.16f %.16f\n', ...
            i,displacements(2*i-1),displacements(2*i),...
            sqrt(displacements(2*i-1)*displacements(2*i-1)+displacements(2*i)*displacements(2*i)));
  end
  % write stresses and strains only if stresses is not empty.
  % since stresses and strains contain these quantities at the subtriangulation
  % level, we only consider the first triangle in the triangulation of each polygon.
  % for this, the helper array triangles_per_polygon is used.
  if ~isempty(stresses)
%     if domainMesh.CPP_remove_holes
%       numel=size(domainMesh.connect(1:domainMesh.npolygons),1);
%     else
      numel=size(domainMesh.connect,1);
%     end    
    fprintf(fid,'ResultName: %s\n','Stresses');   
    fprintf(fid,'NumberOfComponents: %d\n',15);
    fprintf(fid,'ComponentsNames: %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n','e_11','e_22','e_33','e_12','s_11','s_22','s_33','s_12','VM','e_1','e_2','e_3','s_1','s_2','s_3');  
    fprintf(fid,'NumberOfPolygons: %d\n',numel);  
    fprintf(fid,'Polygon e_11 e_22 e_33 e_12 s_11 s_22 s_33 s_12 VM e_1 e_2 e_3 s_1 s_2 s_3\n');    
    k = 1;
    for elem_i=1:numel
      fprintf(fid,'%d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n', ...
              elem_i,strains.e11(k),strains.e22(k),strains.e33(k),...
              strains.e12(k),stresses.s11(k),stresses.s22(k),...
              stresses.s33(k),stresses.s12(k),stresses.vm(k),...
              strains.e1(k),strains.e2(k),strains.e3(k),...
              stresses.s1(k),stresses.s2(k),stresses.s3(k));   
      num_triangles = triangles_per_polygon(elem_i);
      k = k + num_triangles; % points to the first triangle into the next polygon subdivision            
    end
  else
    throw_error('Error in write_solution_CPP_linelast2d.m --> write_solution_CPP_VEM2D_FEM2DT3_linelast2d: stresses container is empty\n');
  end
  fprintf('Check CPP output files in folder: %s\n',...
           config.CPP_output_folder_location);   
  fclose(fid); 
end

function write_polymesh_CPP_VEM2D_linelast2d(domainMesh,config)
  fprintf('\n'); 
  fprintf('Writing polygonal mesh to a text file to be read in Convex Polygon Packing (CPP) program...\n'); 
  % output file
  output_filename=strcat('meshCPP_',config.mesh_filename);
  outfile=[config.CPP_output_folder_location,output_filename];
  fid=fopen(outfile,'w');  
  % write mesh
  num_nodes=size(domainMesh.coords,1);   
%   if domainMesh.CPP_remove_holes
%     numel=size(domainMesh.connect(1:domainMesh.npolygons),1);
%   else
    numel=size(domainMesh.connect,1);
%   end   
  fprintf(fid,'%d %d\n',num_nodes,numel);
  for i = 1:num_nodes
    fprintf(fid,'%.20f %.20f\n',domainMesh.coords(i,1),domainMesh.coords(i,2));
  end
  for i = 1:numel
    connect = domainMesh.connect{i,1};
    nnodes = length(connect);
    fprintf(fid,'%d ',nnodes);
    for j = 1:(nnodes-1)
      fprintf(fid,'%d ',connect(j));
    end
    fprintf(fid,'%d\n',connect(nnodes));
  end

  fprintf('Check CPP output files in folder: %s\n',...
           config.CPP_output_folder_location);   
  fclose(fid); 
end

