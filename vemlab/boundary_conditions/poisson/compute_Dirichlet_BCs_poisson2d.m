function DB_dofs =...
            compute_Dirichlet_BCs_poisson2d(domainMesh,Dirichet_boundary_nodes,...
                                            Dirichet_boundary_dofs,...
                                            Dirichlet_fun_values)
  coords=domainMesh.coords;
  x=coords(Dirichet_boundary_nodes,1); % x-coord of boundary nodes (right boundary)
  y=coords(Dirichet_boundary_nodes,2); % y-coord of boundary nodes (right boundary)    
  %num_nodes=length(x);    
  %DB_dofs.values(1:num_nodes)=Dirichlet_fun_values(x,y);    
  DB_dofs.values=Dirichlet_fun_values(x,y);  
  DB_dofs.indexes=Dirichet_boundary_dofs;      
end