function [K_global,f_global] = assembly(mesh,config,props,body_force_fun_value)
  if strcmp(config.vemlab_module,'LinearElastostatics')
    if strcmp(config.vemlab_method,'VEM2D')
      [K_global,f_global]=vem_assembly_linelast2d(mesh,props,body_force_fun_value);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      [K_global,f_global]=fem_assembly_linelast2d(mesh,config,props,body_force_fun_value);      
    else
      throw_error('In assembly.m: vemlab_method\n');
    end
  else
    throw_error('In assembly.m: vemlab_module\n');
  end
end

