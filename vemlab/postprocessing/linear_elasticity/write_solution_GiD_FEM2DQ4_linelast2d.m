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
% write_solution_GiD_FEM2DQ4_linelast2d(mesh,displacements,stresses,...
%                                       strains,config)
%
% Input
% =====
% mesh          : structure containing the polygonal mesh information
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

function write_solution_GiD_FEM2DQ4_linelast2d(mesh,displacements,stresses,strains,config)

  fprintf('\n'); 
  fprintf('Writing %s solution to a GiD file...\n',config.vemlab_method); 
  output_filename_mesh=strcat(config.mesh_filename,'.flavia.msh');
  outfile_mesh=[config.GiD_output_folder_location,output_filename_mesh]; 
  output_filename_results=strcat(config.mesh_filename,'.flavia.res');
  outfile_results=[config.GiD_output_folder_location,output_filename_results]; 
  
  %% WRITE GiD MESH FILE  
  numel=size(mesh.connect,1);
  nnodes=size(mesh.coords,1);
  fid=fopen(outfile_mesh,'wt'); 
  fprintf(fid,'MESH    dimension 2 ElemType Quadrilateral  Nnode 4\n');
  fprintf(fid,'Coordinates\n');
  for i=1:nnodes
    fprintf(fid,'%d %30.20f %30.20f\n',i,mesh.coords(i,1),mesh.coords(i,2));
  end
  fprintf(fid,'end coordinates\n');
  fprintf(fid,'\n');
  fprintf(fid,'Elements\n');  
  for i=1:numel
    connect=(mesh.connect{i}(:))';
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
  fprintf(fid,'Result  "s11"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "s22"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "s33"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "s12"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "s1"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "s2"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "s3"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "von-Mises"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  %% WRITE STRAINS
  fprintf(fid,'Result  "e11"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "e22"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "e33"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "e12"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "e1"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "e2"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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

  fprintf(fid,'Result  "e3"   "Load Analysis"   %d   Scalar OnGaussPoints "Given gauss points"\n',step);
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
  
  fprintf('Check GiD output files in folder: %s\n',...
           config.GiD_output_folder_location);   
  fclose(fid);
  
end