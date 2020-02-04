function [K_global,f_global] = fem_assembly_poisson2d(domainMesh,config,matProps,source_term_fun_values)
  num_nodes=length(domainMesh.coords(:,1));    
  numel=length(domainMesh.connect);
  K_global=sparse(num_nodes,num_nodes);
  f_global=sparse(num_nodes,1);      
  size_k=length(matProps.k);
  if size_k==numel
    k_is_unique = false;
  elseif size_k==1
    k_is_unique = true;
  end  
  for e=1:numel 
    nodes=domainMesh.connect{e,1};
    verts=domainMesh.coords(nodes,:);
    % conductivity for isotropic material 
    if k_is_unique    
      k=matProps.k; % conductivity is the same for all the elements
    else
      k=matProps.k{e,1}; % conductivity is particular for the current element
    end
    size_k_e=length(k);     
    % assembly element stiffness matrix and element force vector 
    if size_k_e == 1
      K_global(nodes,nodes)=K_global(nodes,nodes)+...
                        fem_stiffness_poisson2d(verts,k,config);
      f_global(nodes)=f_global(nodes)+...
                        fem_source_vector_poisson2d(verts,config,source_term_fun_values); 
    else
      throw_error('In fem_assembly_poisson2d.m: either conductivity was not defined or multiple conductivities assigned to the element... element stiffness not implemented for this condition');
    end    
  end
end

