function DB_dofs =...
            compute_Dirichlet_BCs_linelast2d(mesh,Dirichet_boundary_nodes,...
                                             Dirichet_boundary_dofs,...
                                             Dirichlet_fun_value)
  coords=mesh.coords;
  x=coords(Dirichet_boundary_nodes,1); % x-coord of boundary nodes (right boundary)
  y=coords(Dirichet_boundary_nodes,2); % y-coord of boundary nodes (right boundary)    
  num_nodes=length(x);
  range=1:num_nodes;      
  DB_dofs.values(2*range-1)=Dirichlet_fun_value.ux(x,y);    
  DB_dofs.values(2*range)=Dirichlet_fun_value.uy(x,y);
  DB_dofs.indexes=Dirichet_boundary_dofs;      
end