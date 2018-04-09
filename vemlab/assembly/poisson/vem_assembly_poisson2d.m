function [K_global,f_global] = vem_assembly_poisson2d(domainMesh,matProps,source_term_fun_values)
  num_nodes=length(domainMesh.coords(:,1)); 
  numel=length(domainMesh.connect);  
  K_global=sparse(num_nodes,num_nodes);
  f_global=sparse(num_nodes,1);    
  for e=1:numel 
    nodes=domainMesh.connect{e,1};
    verts=domainMesh.coords(nodes,:);
    % assembly element stiffness matrix and element force vector 
    K_global(nodes,nodes)=K_global(nodes,nodes)+vem_stiffness_poisson2d(verts,matProps);
    f_global(nodes)=f_global(nodes)+vem_source_vector_poisson2d(verts,source_term_fun_values);      
  end
end

