function compute_norms_of_the_error(exact_sol,uh_global,mesh,config,props)
  if strcmp(config.vemlab_module,'LinearElastostatics')
    if strcmp(config.vemlab_method,'VEM2D')
      vem_norms_linelast2d(exact_sol,uh_global,mesh,config,props);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      fem_norms_linelast2d(exact_sol,uh_global,mesh,config,props);      
    else
      throw_error('In compute_norms_of_the_error.m: vemlab_method\n');
    end
  else
    throw_error('In compute_norms_of_the_error.m: vemlab_module\n');
  end
end

