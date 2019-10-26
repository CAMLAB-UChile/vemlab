%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                     VEMLab
%           Source code  : http://camlab.cl/research/software/vemlab/
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% FUNCTION:                PolyMesher2VEMLab
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Read a PolyMesher [1] mesh and write it into a VEMLab mesh format.
%
% Usage
% =====
% PolyMesher2VEMLab(Node,Element,NElem,BoundaryNodes)
%
% Input
% =====
% Node    : PolyMesher array containing the nodal coordinates
% Element : PolyMesher cell array containing the element connectivity 
% NElem   : number of polygonal elements
% BoundaryNodes : Structure containing the boundary nodes
% MeshFile : Location and name of the file for writing the mesh
%
% Output
% ======
%
%-------------------------------------------------------------------------------
% References 
% ==========
% [1] C Talischi, GH Paulino, A Pereira, IFM Menezes, 
%     "PolyMesher: A general-purpose mesh generator for polygonal elements 
%     written in Matlab", Struct Multidisc Optim, 2012,
%     DOI 10.1007/s00158-011-0706-z                                      
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Nov 5, 2018: add plate with a quarter hole domain (by A. Ortiz-Bernardin)
% May 16, 2018: add plate with a hole domain (by A. Ortiz-Bernardin)
% May 15, 2018: add wrench domain (by A. Ortiz-Bernardin)
% Dec. 26, 2017: first realease (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PolyMesher2VEMLab(Node,Element,NElem,BoundaryNodes,MeshFile,...
                           DomainType,BdBox)
  % for now, only rectangular domain is available
  if strcmp(DomainType,'RectangularDomain')
    PolyMesher2VEMLab_rectangular_domain(Node,Element,NElem,BoundaryNodes,...
                                         MeshFile,DomainType);
  elseif strcmp(DomainType,'WrenchDomain')
    PolyMesher2VEMLab_wrench_domain(Node,Element,NElem,BoundaryNodes,...
                                    MeshFile,DomainType,BdBox);  
  elseif strcmp(DomainType,'PlateWithHoleDomain')
    PolyMesher2VEMLab_plate_with_hole_domain(Node,Element,NElem,BoundaryNodes,...
                                             MeshFile,DomainType,BdBox);
  elseif strcmp(DomainType,'PlateWithQuarterHoleDomain')
    PolyMesher2VEMLab_plate_with_quarter_hole_domain(Node,Element,NElem,BoundaryNodes,...
                                                     MeshFile,DomainType,BdBox);                                            
  else
    throw_error('Error in PolyMesher2VEMLab.m --> DomainType\n');
  end
end

function PolyMesher2VEMLab_rectangular_domain(Node,Element,NElem,...
                                              BoundaryNodes,MeshFile,DomainType)
  fprintf('Printing mesh to a VEMLab mesh format...\n'); 
  fid = fopen(MeshFile,'w');
  % print domain type
  fprintf(fid,'# domain type\n');  
  fprintf(fid,'%s\n',DomainType);  
  % print nodal coordinates
  fprintf(fid,'# nodal coordinates: number of nodes followed by the coordinates\n');
  nnode = size(Node,1);
  fprintf(fid,'%d\n',nnode);                                    
  for node_i = 1:nnode
    fprintf(fid,'%.16f %.16f\n', Node(node_i,1), Node(node_i,2));  
  end
  % print element connectivity
  fprintf(fid,'# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel)\n');  
  fprintf(fid,'%d\n',NElem);                                 
  for el = 1:NElem
    NVertex = length(Element{el});
    fprintf(fid,'%d ', NVertex);
    for vertex = 1:(NVertex-1)
      fprintf(fid,'%d ', Element{el}(vertex));
    end
    fprintf(fid,'%d\n', Element{el}(NVertex));
  end
  % print bottom boundary  
  fprintf(fid,'# indices of nodes located on the bottom boundary\n');
  fprintf(fid,'%d ',BoundaryNodes.bottom);  
  fprintf(fid,'\n');
  % print top boundary  
  fprintf(fid,'# indices of nodes located on the top boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.top);
  fprintf(fid,'\n');  
  % print left boundary
  fprintf(fid,'# indices of nodes located on the left boundary\n');
  fprintf(fid,'%d ',BoundaryNodes.left);  
  fprintf(fid,'\n');  
  % print right boundary  
  fprintf(fid,'# indices of nodes located on the right boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.right);    
  fprintf(fid,'\n');
  % print xmin, xmax, ymin, ymax for the rectangular domain
  fprintf(fid,'# xmin, xmax, ymin, ymax of the bounding box\n');  
  xmin=Node(BoundaryNodes.blcorner,1);
  xmax=Node(BoundaryNodes.trcorner,1);
  ymin=Node(BoundaryNodes.blcorner,2);
  ymax=Node(BoundaryNodes.trcorner,2);  
  fprintf(fid,'%.16f %.16f %.16f %.16f\n', xmin, xmax, ymin, ymax);   
  fclose(fid);
end

function PolyMesher2VEMLab_wrench_domain(Node,Element,NElem,BoundaryNodes,...
                                         MeshFile,DomainType,BdBox)
  fprintf('Printing mesh to a VEMLab mesh format...\n'); 
  fid = fopen(MeshFile,'w');
  % print domain type
  fprintf(fid,'# domain type\n');  
  fprintf(fid,'%s\n',DomainType);  
  % print nodal coordinates
  fprintf(fid,'# nodal coordinates: number of nodes followed by the coordinates\n');
  nnode = size(Node,1);
  fprintf(fid,'%d\n',nnode);                                    
  for node_i = 1:nnode
    fprintf(fid,'%.16f %.16f\n', Node(node_i,1), Node(node_i,2));  
  end
  % print element connectivity
  fprintf(fid,'# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel)\n');  
  fprintf(fid,'%d\n',NElem);                                 
  for el = 1:NElem
    NVertex = length(Element{el});
    fprintf(fid,'%d ', NVertex);
    for vertex = 1:(NVertex-1)
      fprintf(fid,'%d ', Element{el}(vertex));
    end
    fprintf(fid,'%d\n', Element{el}(NVertex));
  end
  % print right circle boundary  
  fprintf(fid,'# indices of nodes located on the right circle boundary\n');
  fprintf(fid,'%d ',BoundaryNodes.RightCircle);  
  fprintf(fid,'\n');
  % print left half circle boundary  
  fprintf(fid,'# indices of nodes located on the left half circle boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.LeftHalfCircle);
  fprintf(fid,'\n');  
  % print xmin, xmax, ymin, ymax for the wrench domain
  fprintf(fid,'# xmin, xmax, ymin, ymax of the bounding box\n'); 
  xmin=BdBox(1);
  xmax=BdBox(2);
  ymin=BdBox(3);
  ymax=BdBox(4);  
  fprintf(fid,'%.16f %.16f %.16f %.16f\n', xmin, xmax, ymin, ymax);   
  fclose(fid);
end

function PolyMesher2VEMLab_plate_with_hole_domain(Node,Element,NElem,BoundaryNodes,...
                                                  MeshFile,DomainType,BdBox)
  fprintf('Printing mesh to a VEMLab mesh format...\n'); 
  fid = fopen(MeshFile,'w');
  % print domain type
  fprintf(fid,'# domain type\n');  
  fprintf(fid,'%s\n',DomainType);  
  % print nodal coordinates
  fprintf(fid,'# nodal coordinates: number of nodes followed by the coordinates\n');
  nnode = size(Node,1);
  fprintf(fid,'%d\n',nnode);                                    
  for node_i = 1:nnode
    fprintf(fid,'%.16f %.16f\n', Node(node_i,1), Node(node_i,2));  
  end
  % print element connectivity
  fprintf(fid,'# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel)\n');  
  fprintf(fid,'%d\n',NElem);                                 
  for el = 1:NElem
    NVertex = length(Element{el});
    fprintf(fid,'%d ', NVertex);
    for vertex = 1:(NVertex-1)
      fprintf(fid,'%d ', Element{el}(vertex));
    end
    fprintf(fid,'%d\n', Element{el}(NVertex));
  end
  % print left circle boundary  
  fprintf(fid,'# indices of nodes located on the left boundary\n');
  fprintf(fid,'%d ',BoundaryNodes.Left);  
  fprintf(fid,'\n');
  % print right boundary  
  fprintf(fid,'# indices of nodes located on the right boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.Right);
  fprintf(fid,'\n');   
  % print node located on the middle of the right boundary  
  fprintf(fid,'# index of the node located on the middle of the right boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.MidNodeRightFace);
  fprintf(fid,'\n');    
  % print xmin, xmax, ymin, ymax for the plate with hole domain
  fprintf(fid,'# xmin, xmax, ymin, ymax of the bounding box\n');   
  xmin=BdBox(1);
  xmax=BdBox(2);
  ymin=BdBox(3);
  ymax=BdBox(4);  
  fprintf(fid,'%.16f %.16f %.16f %.16f\n', xmin, xmax, ymin, ymax);   
  fclose(fid);
end

function PolyMesher2VEMLab_plate_with_quarter_hole_domain(Node,Element,NElem,BoundaryNodes,...
                                                          MeshFile,DomainType,BdBox)
  fprintf('Printing mesh to a VEMLab mesh format...\n'); 
  fid = fopen(MeshFile,'w');
  % print domain type
  fprintf(fid,'# domain type\n');  
  fprintf(fid,'%s\n',DomainType);  
  % print nodal coordinates
  fprintf(fid,'# nodal coordinates: number of nodes followed by the coordinates\n');
  nnode = size(Node,1);
  fprintf(fid,'%d\n',nnode);                                    
  for node_i = 1:nnode
    fprintf(fid,'%.16f %.16f\n', Node(node_i,1), Node(node_i,2));  
  end
  % print element connectivity
  fprintf(fid,'# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel)\n');  
  fprintf(fid,'%d\n',NElem);                                 
  for el = 1:NElem
    NVertex = length(Element{el});
    fprintf(fid,'%d ', NVertex);
    for vertex = 1:(NVertex-1)
      fprintf(fid,'%d ', Element{el}(vertex));
    end
    fprintf(fid,'%d\n', Element{el}(NVertex));
  end
  % print left circle boundary  
  fprintf(fid,'# indices of nodes located on the left boundary\n');
  fprintf(fid,'%d ',BoundaryNodes.Left);  
  fprintf(fid,'\n');
  % print right boundary  
  fprintf(fid,'# indices of nodes located on the right boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.Right);
  fprintf(fid,'\n');   
  % print node located on the middle of the right boundary  
  fprintf(fid,'# index of the node located on the middle of the right boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.MidNodeRightFace);
  fprintf(fid,'\n');    
  % print xmin, xmax, ymin, ymax for the plate with hole domain
  fprintf(fid,'# xmin, xmax, ymin, ymax of the bounding box\n');   
  xmin=BdBox(1);
  xmax=BdBox(2);
  ymin=BdBox(3);
  ymax=BdBox(4);  
  fprintf(fid,'%.16f %.16f %.16f %.16f\n', xmin, xmax, ymin, ymax);   
  fclose(fid);
end
