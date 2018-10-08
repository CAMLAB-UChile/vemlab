function write_solution_GiD_poisson2d(domainMesh,scalar_sol,flux,grad,config)
  if strcmp(config.vemlab_method,'VEM2D')||strcmp(config.vemlab_method,'FEM2DT3')
    write_solution_GiD_VEM2D_FEM2DT3_poisson2d(domainMesh,scalar_sol,flux,...
                                                grad,config);
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    write_solution_GiD_FEM2DQ4_poisson2d(domainMesh,scalar_sol,flux,...
                                          grad,config);    
  else
    throw_error('Error in write_solution_GiD_poisson2d.m: vemlab_method')
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:            write_solution_GiD_VEM2D_FEM2DT3_poisson2d 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Write solutions to a GiD [1] output file corresponding to a mesh of 3-node
% triangles. If elements are polygons, they are broken into 3-node triangles.
%
% Usage
% =====
% write_solution_GiD_VEM2D_FEM2DT3_poisson2d(domainMesh,scalar_sol,fluxT3,...
%                                             gradientT3,config)
%
% Input
% =====
% domainMesh    : structure containing the polygonal mesh information
% scalar_sol    : nodal scalar field solution
% fluxT3        : structure storing flux components
% gradientT3    : structure storing gradient components
% config        : structure storing VEMLab configuration options and behavior
%
%
% NOTE: fluxT3 and gradT3 contains entries that are associated with the
%       Gauss points of a subtriangulation of the original mesh composed of
%       polygonal elements. If the original mesh is a three-node triangular mesh
%       then the subtriangulation coincides with the original mesh.
%
%
% Output
% ======
%
%-------------------------------------------------------------------------------
% References 
% ==========
% [1]GiD The personal pre and post processor, https://www.gidhome.com
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Mar. 17, 2018: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_solution_GiD_VEM2D_FEM2DT3_poisson2d(domainMesh,scalar_sol,...
                                                    fluxT3,gradientT3,config)

  fprintf('\n'); 
  fprintf('Writing %s solution to a GiD file...\n',config.vemlab_method); 
  output_filename_mesh=strcat(config.mesh_filename,'.flavia.msh');
  outfile_mesh=[config.GiD_output_folder_location,output_filename_mesh]; 
  output_filename_results=strcat(config.mesh_filename,'.flavia.res');
  outfile_results=[config.GiD_output_folder_location,output_filename_results]; 
  
  %% WRITE GiD MESH FILE  
  num_polygons=size(domainMesh.connect,1);
  nnodes=size(domainMesh.coords,1);
  fid=fopen(outfile_mesh,'wt'); 
  fprintf(fid,'MESH    dimension 2 ElemType Triangle  Nnode 3\n');
  fprintf(fid,'Coordinates\n');
  for i=1:nnodes
    fprintf(fid,'%d %30.20f %30.20f\n',i,domainMesh.coords(i,1),domainMesh.coords(i,2));
  end
  fprintf(fid,'end coordinates\n');
  fprintf(fid,'\n');
  fprintf(fid,'Elements\n');
  k=0;
  for i=1:num_polygons
    % make the subtriangulation of the polygonal mesh that will be used with
    % the entries stored in fluxT3 and gradientT3 --- see note above.    
    if strcmp(config.vemlab_method,'VEM2D')
      connect=triangulate_polygon(domainMesh,i);
    elseif strcmp(config.vemlab_method,'FEM2DT3')
      node_indices=domainMesh.connect(i,:);
      connect=(node_indices{1}(:))';
    else
      throw_error('Error in write_solution_GiD_poisson2d.m --> write_solution_GiD_VEM2D_FEM2DT3_poisson2d: vemlab_method\n');
    end
    for tr_i=1:size(connect,1)
      k=k+1;   
      % write mesh of triangles
      fprintf(fid,'%d %d %d %d %d\n',k,connect(tr_i,1),connect(tr_i,2),connect(tr_i,3),1);  
    end
  end
  num_triangles=k;
  fprintf(fid,'end elements\n');
  fclose(fid);

  %% WRITE INFO FOR PLOTTING ON GAUSS POINTS
  npoints=1;
  fid=fopen(outfile_results,'wt'); 
  fprintf(fid,'Gid Post Results File 1.0\n');
  fprintf(fid,'GaussPoints "Given gauss points" ElemType Triangle\n');
  fprintf(fid,'Number of Gauss Points: %d\n',npoints);
  fprintf(fid,'Natural Coordinates: Given\n');
  tcoord=[1.0/3.0 1.0/3.0 1.0/3.0];  
  for i=1:npoints
    fprintf(fid,'%f %f\n',tcoord(i,1),tcoord(i,2));
  end      
  fprintf(fid,'End gausspoints\n');   
  
  %% WRITE SCALAR SOLUTION
  step=1;
  fprintf(fid,'Result  "Scalar field//Scalar Field"   "Load Analysis"   %d   Scalar OnNodes\n',step);
  fprintf(fid,'Values\n');
  for i=1:nnodes
    fprintf(fid,'%d %30.20f\n',i,scalar_sol(i));
  end
  fprintf(fid,'End Values\n');  
  
  %% WRITE FLUX
  
  % write fluxes only if fluxT3 is not empty
  if ~isempty(fluxT3)    
    fluxNorm=sqrt((fluxT3.qx).*(fluxT3.qx)+(fluxT3.qy).*(fluxT3.qy));      
    
    fprintf(fid,'Result  "Fluxes//Flux-x"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',fluxT3.qx(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Fluxes//Flux-y"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',fluxT3.qy(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n'); 
    
    fprintf(fid,'Result  "Fluxes//||Flux||"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',fluxNorm(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');     
  end

  %% WRITE GRADIENT
  
  % write gradient only if gradientT3 is not empty
  if ~isempty(gradientT3)   
    gradNorm=sqrt((gradientT3.dx).*(gradientT3.dx)+(gradientT3.dy).*(gradientT3.dy));  
    
    fprintf(fid,'Result  "Gradient//Grad-x"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',gradientT3.dx(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Gradient//Grad-y"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',gradientT3.dy(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n'); 
    
    fprintf(fid,'Result  "Gradient//||Grad||"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',gradNorm(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');      
  end
  fprintf('Check GiD output files in folder: %s\n',...
           config.GiD_output_folder_location);   
  fclose(fid);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:             write_solution_GiD_FEM2DQ4_poisson2d 
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Write solutions to a GiD [1] output file corresponding to a mesh of 4-node
% quadrilaterals. For plotting of stresses and strains, 2x2 quadrature rule is
% used.
%
% Usage
% =====
% write_solution_GiD_FEM2DQ4_poisson2d(domainMesh,scalar_sol,flux,grad,config)
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
% [1]GiD The personal pre and post processor, https://www.gidhome.com
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Mar. 17, 2018: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_solution_GiD_FEM2DQ4_poisson2d(domainMesh,scalar_sol,flux,grad,config)

  fprintf('\n'); 
  fprintf('Writing %s solution to a GiD file...\n',config.vemlab_method); 
  output_filename_mesh=strcat(config.mesh_filename,'.flavia.msh');
  outfile_mesh=[config.GiD_output_folder_location,output_filename_mesh]; 
  output_filename_results=strcat(config.mesh_filename,'.flavia.res');
  outfile_results=[config.GiD_output_folder_location,output_filename_results]; 
  
  %% WRITE GiD MESH FILE  
  numel=size(domainMesh.connect,1);
  nnodes=size(domainMesh.coords,1);
  fid=fopen(outfile_mesh,'wt'); 
  fprintf(fid,'MESH    dimension 2 ElemType Quadrilateral  Nnode 4\n');
  fprintf(fid,'Coordinates\n');
  for i=1:nnodes
    fprintf(fid,'%d %30.20f %30.20f\n',i,domainMesh.coords(i,1),domainMesh.coords(i,2));
  end
  fprintf(fid,'end coordinates\n');
  fprintf(fid,'\n');
  fprintf(fid,'Elements\n');  
  for i=1:numel
    connect=(domainMesh.connect{i}(:))';
    fprintf(fid,'%d %d %d %d %d %d\n',i,connect(1),connect(2),connect(3),connect(4),1);
  end
  fprintf(fid,'end elements\n');  
  fclose(fid);

  %% WRITE INFO FOR PLOTTING ON GAUSS POINTS
  npoints=4; % 2x2 quadrature rule is used
  fid=fopen(outfile_results,'wt'); 
  fprintf(fid,'Gid Post Results File 1.0\n');
  fprintf(fid,'GaussPoints "Given gauss points" ElemType Quadrilateral\n');
  fprintf(fid,'Number of Gauss Points: %d\n',npoints);
  fprintf(fid,'Natural Coordinates: Given\n');
  xi(1,1) = -0.5773502692;  xi(2,1) = xi(1,1);
  xi(1,2) = -xi(1,1);       xi(2,2) = xi(1,1);
  xi(1,3) = xi(1,1);        xi(2,3) = -xi(1,1);
  xi(1,4) = -xi(1,1);       xi(2,4) = -xi(1,1);
  for i=1:npoints
    fprintf(fid,'%f %f\n',xi(1,i),xi(2,i));
  end        
  fprintf(fid,'End gausspoints\n');      
  
  %% WRITE SCALAR SOLUTION
  step=1;
  fprintf(fid,'Result  "Scalar field//Scalar Field"   "Load Analysis"   %d   Scalar OnNodes\n',step);
  fprintf(fid,'Values\n');
  for i=1:nnodes
    fprintf(fid,'%d %30.20f\n',i,scalar_sol(i));
  end
  fprintf(fid,'End Values\n');  
  
  %% WRITE FLUX
  
  % write fluxes only if flux is not empty
  if ~isempty(flux)      
    fluxNorm=sqrt((flux.qx).*(flux.qx)+(flux.qy).*(flux.qy));  
    
    fprintf(fid,'Result  "Fluxes//Flux-x"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',flux.qx(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Fluxes//Flux-y"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',flux.qy(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');  
    
    fprintf(fid,'Result  "Fluxes//||Flux||"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',fluxNorm(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');     
  end

  %% WRITE GRADIENT
  
  % write gradients only if grad is not empty
  if ~isempty(grad)  
    gradNorm=sqrt((grad.dx).*(grad.dx)+(grad.dy).*(grad.dy)); 
    
    fprintf(fid,'Result  "Gradient//Grad-x"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',grad.dx(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Gradient//Grad-y"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',grad.dy(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');  
    
    fprintf(fid,'Result  "Gradient//||Grad||"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',gradNorm(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');     
  end
  fprintf('Check GiD output files in folder: %s\n',...
           config.GiD_output_folder_location);   
  fclose(fid);
  
end

