function DB_dofs = compute_Dirichlet_BCs(mesh,config,Dirichet_boundary_nodes,...
                                         Dirichet_boundary_dofs,Dirichlet_fun_value)
  if strcmp(config.vemlab_module,'LinearElastostatics')
    if strcmp(config.vemlab_method,'VEM2D')
      DB_dofs=compute_Dirichlet_BCs_linelast2d(mesh,Dirichet_boundary_nodes,...
                                               Dirichet_boundary_dofs,...
                                               Dirichlet_fun_value);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      DB_dofs=compute_Dirichlet_BCs_linelast2d(mesh,Dirichet_boundary_nodes,...
                                               Dirichet_boundary_dofs,...
                                               Dirichlet_fun_value);      
    else
      throw_error('In compute_Dirichlet_BCs.m: vemlab_method\n');
    end
  else
    throw_error('In compute_Dirichlet_BCs.m: vemlab_module\n');
  end
end

