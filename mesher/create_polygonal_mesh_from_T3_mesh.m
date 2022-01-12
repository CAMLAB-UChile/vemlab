function create_polygonal_mesh_from_T3_mesh
  close all;
  clc;
  
  %% ADD ALL VEMLAB FOLDERS TO THE PATH
  opsystem=computer;
  is_Windows = strcmp(opsystem,'PCWIN') || strcmp(opsystem,'PCWIN64');
  is_Linux = strcmp(opsystem,'GLNX86') || strcmp(opsystem,'GLNXA64');   
  if is_Windows
    cd ..\; vemlab_root_dir=setpath;
  elseif is_Linux
    cd ../; vemlab_root_dir=setpath;
  end   

  %%
  %%%%%%%%%%%%%%%%%%           USER AREA         %%%%%%%%%%%%%%%%%%%%%%%%%
  
  input_T3_mesh_filename='GiD_T3_thick_walled_cylinder_40x80.txt';
  output_polygonal_mesh_filename='thick_walled_cylinder_poly_elems.txt';  
  
  % WARNING: DON'T CHANGE THE FOLLOWING TWO LINES IF YOU DON'T KNOW WHAT YOU ARE DOING!
  config=config_vemlab_mesher(opsystem,vemlab_root_dir,input_T3_mesh_filename); % configure mesher
  mesh_file=[config.GiD_T3_mesh_folder_location,input_T3_mesh_filename];
  
  % READ T3 MESH
  domainMesh=read_T3_mesh(mesh_file);
  
  % CREATE POLYGONAL MESH FROM T3 MESH
  Element=domainMesh.connect;
  Nodes=domainMesh.coords;  
  [Element,Node,Dirichlet_nodes,Neumann_nodes]=...
     T3toPolyMesh(Element,Nodes,domainMesh.boundary_edges,...
                  domainMesh.corner_nodes,domainMesh.num_dirichlet_entities,...
                  domainMesh.num_neumann_entities);
  
%   % distort the mesh (comment if not required)
%   b=3.4;
%   a=0;
%   Node(:,1)=Node(:,1)+0.3*Node(:,1).*sin(b*Node(:,1)+a).*sin(b*Node(:,2)+a).*(BdBox(2)-Node(:,1))./BdBox(2);
%   Node(:,2)=Node(:,2)+0.3*Node(:,2).*sin(b*Node(:,2)+a).*sin(b*Node(:,1)+a).*(BdBox(4)-Node(:,2))./BdBox(4);

  %Plot mesh to a VEMLab mesh format
  mesh_file=[config.mesh_folder_location,output_polygonal_mesh_filename];
  write_vemlab_meshfile(Element,Node,Dirichlet_nodes,Neumann_nodes,mesh_file,domainMesh.domain_type,domainMesh.BdBox);      
 
  % plot mesh
  config.plot_mesh = 'yes'; % 'yes' or 'no'
  config.plot_node_numbers = 'yes'; % 'yes' or 'no'
  config.plot_element_numbers = 'no'; % 'yes' or 'no'  
  config.plot_axis = 'yes'; % 'yes' or 'no'
  %%%%%%%%%%%%%%%%%%         END USER AREA         %%%%%%%%%%%%%%%%%%%%%%%  
  plot_mesh(Element,Node,config);  
  
end

function write_vemlab_meshfile(Element,Node,DirichletNodes,NeumannNodes,MeshFile,...
                               DomainType,BdBox)
  % for now, only rectangular domain is available
  if strcmp(DomainType,'CustomBoundaryDataWithTags')
    write_custom_vemlab_meshfile_from_T3(Element,Node,DirichletNodes,NeumannNodes,...
                                        MeshFile,DomainType,BdBox);                                          
  else
    throw_error('Error in create_polygonal_mesh_from_T3_mesh.m --> wrte_vemlab_meshfile.m --> DomainType\n');
  end
end

function write_custom_vemlab_meshfile_from_T3(Element,Node,DirichletNodes,NeumannNodes,...
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
  NElem = length(Element);
  fprintf(fid,'# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 node 3 node 4)\n');  
  fprintf(fid,'%d\n',NElem);                                 
  for el = 1:NElem
    NVertex = length(Element{el});
    fprintf(fid,'%d ', NVertex);
    for vertex = 1:(NVertex-1)
      fprintf(fid,'%d ', Element{el}(vertex));
    end
    fprintf(fid,'%d\n', Element{el}(NVertex));
  end
  % Dirichlet boundary nodes
  fprintf(fid,'# indices of nodes located on the Dirichlet boundary: first line, Dirichlet nodes, second line, tags associated with the Dirichlet nodes\n');
  if isempty(DirichletNodes(:,1))
    fprintf(fid,'%d ',0);    
    fprintf(fid,'\n');
    fprintf(fid,'%d ',0);
  else 
    fprintf(fid,'%d ',DirichletNodes(:,1));
    fprintf(fid,'\n');    
    fprintf(fid,'%d ',DirichletNodes(:,2));    
  end
  fprintf(fid,'\n');
  % Neumann boundary nodes
  fprintf(fid,'# indices of nodes located on the Neumann boundary: first line, Neumann nodes, second line, tags associated with the Neumann nodes. ***WARNING: MUST BE ORDERED BY HAND IF NODES ARE NOT GIVEN IN CONSECUTIVE ORDER***\n'); 
  if isempty(NeumannNodes(:,1))
    fprintf(fid,'%d ',0);
    fprintf(fid,'\n');    
    fprintf(fid,'%d ',0);    
  else
    fprintf(fid,'%d ',NeumannNodes(:,1));
    fprintf(fid,'\n');    
    fprintf(fid,'%d ',NeumannNodes(:,2));    
  end
  fprintf(fid,'\n');  
  % print xmin, xmax, ymin, ymax for the rectangular domain
  fprintf(fid,'# xmin, xmax, ymin, ymax of the bounding box\n');  
  fprintf(fid,'%.16f %.16f %.16f %.16f\n',BdBox(1),BdBox(2),BdBox(3),BdBox(4));   
  fclose(fid);
end

function plot_mesh(connect,coords,config)

  if strcmp(config.plot_mesh,'yes')
    nelems=length(connect);
    polygons=cell(nelems,1);
    for e=1:nelems
      polygons{e,1}=connect{e,1};
    end 
    figure; 
    maxNumVertices = max(cellfun(@numel,polygons));
    padFunc = @(vertList) [vertList' NaN(1,maxNumVertices-numel(vertList))];
    elements = cellfun(padFunc,polygons,'UniformOutput',false);
    elements = vertcat(elements{:}); 
    patch('Faces',elements,'Vertices',coords,'FaceColor','w'); 
    axis equal; 
    if strcmp(config.plot_axis,'yes')
     axis on;
    else
      axis off;
    end

    hold on;

    % set the font size for the node and element numbering
    if nelems <= 50
      fsize = 14; % 10
    else
      fsize = 10; % 6
    end
    
    dx = 0; dy = 0;
    xlim([min(coords(:, 1)) - dx, max(coords(:, 1)) + dx]);
    ylim([min(coords(:, 2)) - dy, max(coords(:, 2)) + dy]);
    xlabel('$x$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');
    ylabel('$y$','FontWeight','bold','FontSize',20,'FontName',...
           'Times New Roman','Interpreter','latex');    

    % node numbering
    numnod=size(coords,1);
    mesh_size=min_edge_size(coords,connect,nelems);
    if strcmp(config.plot_node_numbers,'yes')
      offset_x = 0.06*mesh_size; offset_y = 0.08*mesh_size;
      for i = 1:numnod
        idstring = sprintf('%d',i);
        text(coords(i,1)+offset_x, coords(i,2)+offset_y, idstring,'color',[1,0,0],'fontsize',fsize);
      end
    end

    % element numbering for QUAD4
    if strcmp(config.plot_element_numbers,'yes')
      for e = 1:length(connect(:,1))
        xc = mean(coords(connect{e},1)); % element number is plotted in the
        yc = mean(coords(connect{e},2)); % center of the rectangle
        idstring = sprintf('%d',e);
        text(xc, yc, idstring, 'color',[0,0,0],'fontsize',fsize);
      end
    end  
  end
  
end

function h_min = min_edge_size(Node,Element,NElem)
  h_min=1e1000;
  for e=1:NElem
    NVertex = length(Element{e});
    nodes = Element{e}(1:NVertex);
    verts=Node(nodes,:); 

    mysize=size(verts,1);
    mydist=zeros(mysize*mysize,1);

    % This procedure also compute the distance between a node with itself (0); however, it
    % doesn't matter since at the end the maximum distance is selected
    r = 1;
    for node_i = 1:mysize
      for node_j = 1:mysize
        % shifts coordinates with respect to "node i" in the connectivity
        shifted_coords= verts-repmat(verts(node_i,:),mysize,1);
        % compute the distances wrt to "node_i" in the element connectivity
        for node = 1:(mysize-1)
          mydist(r) = norm(shifted_coords(node,:));
          r = r + 1;
        end
      end  
    end

    % minimum distance from the first node to any node in the element connectivity 
    h_size=min(mydist);

    if h_size<h_min
      h_min=h_size;
    end   
  end
end



