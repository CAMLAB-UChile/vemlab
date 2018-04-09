function [K_global,f_global] = vem_assembly_linelast2d(domainMesh,matProps,body_force_fun_values)
  num_nodes=length(domainMesh.coords(:,1)); 
  numel=length(domainMesh.connect);  
  K_global=sparse(2*num_nodes,2*num_nodes);
  f_global=sparse(2*num_nodes,1);    
  for e=1:numel 
    nodes=domainMesh.connect{e,1};
    verts=domainMesh.coords(nodes,:);
    % assembly element stiffness matrix and element force vector 
    ind=zeros(2*length(nodes),1);      
    range=1:length(nodes);
    ind(2*range-1)=2*nodes-1;
    ind(2*range)=2*nodes;
    K_global(ind,ind)=K_global(ind,ind)+vem_stiffness_linelast2d(verts,matProps);
    f_global(ind)=f_global(ind)+vem_body_force_linelast2d(verts,body_force_fun_values);      
  end
end

