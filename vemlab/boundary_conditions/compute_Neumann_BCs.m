function Neumann_BCs = compute_Neumann_BCs(mesh,config,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_value)
  if strcmp(config.vemlab_module,'LinearElastostatics')
    if strcmp(config.vemlab_method,'VEM2D')
      Neumann_BCs=vem_compute_Neumann_BCs_linelast2d(mesh,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_value);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      Neumann_BCs=fem_compute_Neumann_BCs_linelast2d(mesh,Neumann_boundary_nodes,Neumann_boundary_dofs,Neumann_fun_value);      
    else
      throw_error('In compute_Neumann_BCs.m: vemlab_method\n');
    end
  else
    throw_error('In compute_Neumann_BCs.m: vemlab_module\n');
  end
end

