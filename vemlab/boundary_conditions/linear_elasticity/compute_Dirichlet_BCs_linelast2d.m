function DB_dofs =...
            compute_Dirichlet_BCs_linelast2d(domainMesh,Dirichlet_boundary_nodes,...
                                             Dirichet_boundary_dofs,...
                                             Dirichlet_fun_values,matProps)
  coords=domainMesh.coords;
  x=coords(Dirichlet_boundary_nodes,1); % x-coord of Dirichlet boundary nodes
  y=coords(Dirichlet_boundary_nodes,2); % y-coord of Dirichlet boundary nodes   
  num_DB_nodes=length(x);
  range=1:num_DB_nodes;     
  DB_dofs_aux.values(2*range-1,:)=Dirichlet_fun_values.ux(x,y,matProps);    
  DB_dofs_aux.values(2*range,:)=Dirichlet_fun_values.uy(x,y,matProps);
  DB_dofs_aux.indexes=Dirichet_boundary_dofs;    
  
  % disregard dofs and values of those nodes that are free
  ind=find(DB_dofs_aux.values(:,2)==0);
  DB_dofs.values=DB_dofs_aux.values(ind,1);
  DB_dofs.indexes=DB_dofs_aux.indexes(ind);
end