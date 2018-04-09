function [K_global,f_global] = fem_assembly_poisson2d(domainMesh,config,matProps,source_term_fun_values)
  num_nodes=length(domainMesh.coords(:,1));    
  numel=length(domainMesh.connect);
  K_global=sparse(num_nodes,num_nodes);
  f_global=sparse(num_nodes,1);      
  for e=1:numel 
    nodes=domainMesh.connect{e,1};
    verts=domainMesh.coords(nodes,:);
    % assembly element stiffness matrix and element force vector 
    K_global(nodes,nodes)=K_global(nodes,nodes)+...
                      fem_stiffness_poisson2d(verts,matProps,config);
    f_global(nodes)=f_global(nodes)+...
                      fem_source_vector_poisson2d(verts,config,source_term_fun_values);      
  end
end

