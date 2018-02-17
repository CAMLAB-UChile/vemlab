function [K_global,f_global] = fem_assembly_linelast2d(mesh,config,matProps,body_force_fun_value)
  num_nodes=length(mesh.coords(:,1));    
  numel=length(mesh.connect);
  K_global=sparse(2*num_nodes,2*num_nodes);
  f_global=sparse(2*num_nodes,1);      
  for e=1:numel 
    nodes=mesh.connect{e,1};
    verts=mesh.coords(nodes,:);
    % assembly element stiffness matrix and element force vector 
    ind=zeros(2*length(nodes),1);      
    range=1:length(nodes);
    ind(2*range-1)=2*nodes-1;
    ind(2*range)=2*nodes;
    K_global(ind,ind)=K_global(ind,ind)+...
                      fem_stiffness_linelast2d(verts,matProps,config);
    f_global(ind)=f_global(ind)+...
                      fem_body_force_linelast2d(verts,config,body_force_fun_value);      
  end
end

