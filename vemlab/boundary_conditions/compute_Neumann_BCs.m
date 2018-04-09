function Neumann_BCs = compute_Neumann_BCs(domainMesh,config,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_values)
  if strcmp(config.vemlab_module,'LinearElastostatics')
    if strcmp(config.vemlab_method,'VEM2D')
      Neumann_BCs=vem_compute_Neumann_BCs_linelast2d(domainMesh,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_values);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      Neumann_BCs=fem_compute_Neumann_BCs_linelast2d(domainMesh,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_values);      
    else
      throw_error('In compute_Neumann_BCs.m: vemlab_method\n');
    end
  elseif strcmp(config.vemlab_module,'Poisson')
    if strcmp(config.vemlab_method,'VEM2D')
      Neumann_BCs=vem_compute_Neumann_BCs_poisson2d(domainMesh,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_values);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      Neumann_BCs=fem_compute_Neumann_BCs_poisson2d(domainMesh,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_values);      
    else
      throw_error('In compute_Neumann_BCs.m: vemlab_method\n');
    end    
  else
    throw_error('In compute_Neumann_BCs.m: vemlab_module\n');
  end
end

