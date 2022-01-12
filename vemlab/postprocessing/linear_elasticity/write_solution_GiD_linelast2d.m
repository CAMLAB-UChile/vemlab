function write_solution_GiD_linelast2d(domainMesh,displacements,stress,strain,config)

  if strcmp(config.vemlab_method,'VEM2D')||strcmp(config.vemlab_method,'FEM2DT3')
    write_solution_GiD_VEM2D_FEM2DT3_linelast2d(domainMesh,displacements,stress,...
                                                strain,config);
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    write_solution_GiD_FEM2DQ4_linelast2d(domainMesh,displacements,stress,...
                                          strain,config);    
  else
    throw_error('Error in write_solution_GiD_linelast2d.m: vemlab_method')
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:            write_solution_GiD_VEM2D_FEM2DT3_linelast 
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
% write_solution_GiD_VEM2D_FEM2DT3_linelast2d(domainMesh,displacements,stressesT3,...
%                                             strainsT3,config)
%
% Input
% =====
% domainMesh    : structure containing the polygonal mesh information
% displacements : nodal displacement solution
% stressesT3    : structure storing stress tensor components and vMises stress
% strainsT3     : structure storing strain tensor components
% config        : structure storing VEMLab configuration options and behavior
%
% NOTE: stressesT3 and strainsT3 contains entries that are associated with the
%       Gauss points of a subtriangulation of the original mesh composed of
%       polygonal elements. If the original mesh is a three-node triangular mesh
%       then the subtriangulation coincides with the original mesh.
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
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_solution_GiD_VEM2D_FEM2DT3_linelast2d(domainMesh,displacements,...
                                                     stresses,strains,config)

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
    % the entries stored in stressesT3 and strainsT3 --- see note above.
    if strcmp(config.vemlab_method,'VEM2D')||strcmp(config.vemlab_method,'NIVEM2D')
      connect=triangulate_polygon(domainMesh,i);
    elseif strcmp(config.vemlab_method,'FEM2DT3')
      node_indices=domainMesh.connect(i,:);
      connect=(node_indices{1}(:))'; 
    else
      throw_error('Error in write_solution_GiD_linelast2d.m --> write_solution_GiD_VEM2D_FEM2DT3_linelast2d: vemlab_method\n');
    end
    for tr_i=1:size(connect,1)
      k=k+1;   
      % write mesh of triangles
      fprintf(fid,'%d %d %d %d %d\n',k,connect(tr_i,1),connect(tr_i,2),connect(tr_i,3),1);
      % AOB: (19-AUG-2021) assign the element stresses and strains to the
      % element's subtriangulation since stresses and strains contain one
      % constant value per polygon 
      stressesT3.s11(k)=stresses.s11(i);
      stressesT3.s22(k)=stresses.s22(i);  
      stressesT3.s33(k)=stresses.s33(i); 
      stressesT3.s12(k)=stresses.s12(i);  
      stressesT3.s1(k)=stresses.s1(i);       
      stressesT3.s2(k)=stresses.s2(i);  
      stressesT3.s3(k)=stresses.s3(i);   
      stressesT3.vm(k)=stresses.vm(i);    
      stressesT3.p(k)=stresses.p(i);
      strainsT3.e11(k)=strains.e11(i);    
      strainsT3.e22(k)=strains.e22(i);  
      strainsT3.e33(k)=strains.e33(i);         
      strainsT3.e12(k)=strains.e12(i);
      strainsT3.e1(k)=strains.e1(i);   
      strainsT3.e2(k)=strains.e2(i);  
      strainsT3.e3(k)=strains.e3(i);     
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
  
  %% WRITE DISPLACEMENTS
  step=1;
  fprintf(fid,'Result  "Displacements"   "Load Analysis"   %d   Vector OnNodes\n',step);
  fprintf(fid,'ComponentNames  "X-Disp"   "Y-Disp"\n');
  fprintf(fid,'Values\n');
  for i=1:nnodes
    fprintf(fid,'%d %30.20f %30.20f\n',i,displacements(2*i-1),displacements(2*i));
  end
  fprintf(fid,'End Values\n');  
  
  %% WRITE STRESSES
  
  % write stresses only if stressesT3 is not empty
  if ~isempty(stressesT3)  
    fprintf(fid,'Result  "Stresses//s11"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',stressesT3.s11(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Stresses//s22"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stressesT3.s22(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n'); 

    fprintf(fid,'Result  "Stresses//s33"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stressesT3.s33(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');       

    fprintf(fid,'Result  "Stresses//s12"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stressesT3.s12(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n'); 

    fprintf(fid,'Result  "Principal stresses//s1"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',stressesT3.s1(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Principal stresses//s2"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stressesT3.s2(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');  

    fprintf(fid,'Result  "Principal stresses//s3"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stressesT3.s3(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');       

    fprintf(fid,'Result  "Stresses//von-Mises"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stressesT3.vm(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');  
    
    fprintf(fid,'Result  "Stresses//Pressure"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stressesT3.p(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');     
  end

  %% WRITE STRAINS
  
  % write stresses only if strainsT3 is not empty
  if ~isempty(strainsT3)   
    fprintf(fid,'Result  "Strains//e11"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',strainsT3.e11(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Strains//e22"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',strainsT3.e22(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n'); 

    fprintf(fid,'Result  "Strains//e33"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',strainsT3.e33(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');       

    fprintf(fid,'Result  "Strains//e12"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',strainsT3.e12(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');    

    fprintf(fid,'Result  "Principal strains//e1"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',strainsT3.e1(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Principal strains//e2"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',strainsT3.e2(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n'); 

    fprintf(fid,'Result  "Principal strains//e3"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:num_triangles
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',strainsT3.e3(kk));
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
% FUNCTION:             write_solution_GiD_FEM2DQ4_linelast 
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
% write_solution_GiD_FEM2DQ4_linelast2d(domainMesh,displacements,stresses,...
%                                       strains,config)
%
% Input
% =====
% domainMesh    : structure containing the polygonal mesh information
% displacements : nodal displacement solution
% stresses      : structure storing stress tensor components and vMises stress
% strains       : structure storing strain tensor components
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
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_solution_GiD_FEM2DQ4_linelast2d(domainMesh,displacements,stresses,strains,config)

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
  
  %% WRITE DISPLACEMENTS
  step=1;
  fprintf(fid,'Result  "Displacements"   "Load Analysis"   %d   Vector OnNodes\n',step);
  fprintf(fid,'ComponentNames  "X-Disp"   "Y-Disp"\n');
  fprintf(fid,'Values\n');
  for i=1:nnodes
    fprintf(fid,'%d %30.20f %30.20f\n',i,displacements(2*i-1),displacements(2*i));
  end
  fprintf(fid,'End Values\n');  
  
  %% WRITE STRESSES
  
  % write stresses only if stresses is not empty
  if ~isempty(stresses)   
    fprintf(fid,'Result  "Stresses//s11"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',stresses.s11(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Stresses//s22"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stresses.s22(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n'); 

    fprintf(fid,'Result  "Stresses//s33"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stresses.s33(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');       

    fprintf(fid,'Result  "Stresses//s12"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stresses.s12(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n'); 

    fprintf(fid,'Result  "Principal stresses//s1"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',stresses.s1(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Principal stresses//s2"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stresses.s2(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');  

    fprintf(fid,'Result  "Principal stresses//s3"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stresses.s3(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');       

    fprintf(fid,'Result  "Stresses//von-Mises"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stresses.vm(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');
    
    fprintf(fid,'Result  "Stresses//Pressure"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',stresses.p(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');    
  end

  %% WRITE STRAINS
  
  % write strains only if strains is not empty
  if ~isempty(strains)   
    fprintf(fid,'Result  "Strains//e11"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',strains.e11(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Strains//e22"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',strains.e22(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n'); 

    fprintf(fid,'Result  "Strains//e33"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',strains.e33(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');       

    fprintf(fid,'Result  "Strains//e12"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',strains.e12(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');    

    fprintf(fid,'Result  "Principal strains//e1"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number   
      for j=1:npoints
        fprintf(fid,'%f\n',strains.e1(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');

    fprintf(fid,'Result  "Principal strains//e2"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',strains.e2(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n'); 

    fprintf(fid,'Result  "Principal strains//e3"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
    fprintf(fid,'Values\n');
    kk=1;
    for i=1:numel
      fprintf(fid,'%d ',i); % element number
      for j=1:npoints
        fprintf(fid,'%f\n',strains.e3(kk));
        kk=kk+1;
      end
    end
    fprintf(fid,'End Values\n');
  end
  
  fprintf('Check GiD output files in folder: %s\n',...
           config.GiD_output_folder_location);   
  fclose(fid);
  
end
