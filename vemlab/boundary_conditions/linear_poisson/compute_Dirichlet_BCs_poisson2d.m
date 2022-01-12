%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: compute_Dirichlet_BCs_poisson2d
% Purpose: compute Dirichlet BCs as a function of the nodal coordinates (x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DB_dofs =...
            compute_Dirichlet_BCs_poisson2d(domainMesh,Dirichet_boundary_nodes,...
                                            Dirichet_boundary_dofs,...
                                            Dirichlet_fun_values,matProps)                                                                             
  coords=domainMesh.coords;
  x=coords(Dirichet_boundary_nodes,1); % x-coord of boundary nodes (right boundary)
  y=coords(Dirichet_boundary_nodes,2); % y-coord of boundary nodes (right boundary)    
       
  DB_dofs_aux.values=Dirichlet_fun_values(x,y,matProps);  
  DB_dofs_aux.indexes=Dirichet_boundary_dofs; 
  
  % disregard dofs and values of those nodes that are free
  ind=find(DB_dofs_aux.values(:,2)==0);
  DB_dofs.values=DB_dofs_aux.values(ind,1);
  DB_dofs.indexes=DB_dofs_aux.indexes(ind);  
end