function [K_global,f_global] = assembly(domainMesh,config,props,source_or_body_force_fun_values)
  if strcmp(config.vemlab_module,'LinearElastostatics')
    if strcmp(config.vemlab_method,'VEM2D') && strcmp(config.vemlab_solver,'dense')
      [K_global,f_global]=vem_assembly_linelast2d(domainMesh,props,source_or_body_force_fun_values);
    elseif strcmp(config.vemlab_method,'VEM2D') && strcmp(config.vemlab_solver,'sparse')
      [K_global,f_global]=vem_sparse_assembly_linelast2d(domainMesh,props,source_or_body_force_fun_values);      
    elseif (strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')) && strcmp(config.vemlab_solver,'dense')
      [K_global,f_global]=fem_assembly_linelast2d(domainMesh,config,props,source_or_body_force_fun_values);   
    elseif (strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')) && strcmp(config.vemlab_solver,'sparse')
      [K_global,f_global]=fem_sparse_assembly_linelast2d(domainMesh,config,props,source_or_body_force_fun_values);      
    else
      throw_error('In assembly.m: vemlab_method\n');
    end
  elseif strcmp(config.vemlab_module,'Poisson')
    if strcmp(config.vemlab_method,'VEM2D') && strcmp(config.vemlab_solver,'dense')
      [K_global,f_global]=vem_assembly_poisson2d(domainMesh,props,source_or_body_force_fun_values);
    elseif strcmp(config.vemlab_method,'VEM2D') && strcmp(config.vemlab_solver,'sparse')
      [K_global,f_global]=vem_sparse_assembly_poisson2d(domainMesh,props,source_or_body_force_fun_values);       
    elseif (strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')) && strcmp(config.vemlab_solver,'dense')
      [K_global,f_global]=fem_assembly_poisson2d(domainMesh,config,props,source_or_body_force_fun_values); 
    elseif (strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')) && strcmp(config.vemlab_solver,'sparse')
      [K_global,f_global]=fem_sparse_assembly_poisson2d(domainMesh,config,props,source_or_body_force_fun_values);      
    else
      throw_error('In assembly.m: vemlab_method\n');
    end    
  else
    throw_error('In assembly.m: vemlab_module\n');
  end
end

