function [stresses,strains]=compute_stresses_and_strains(u_nodal_sol,mesh,props,config)
  if strcmp(config.vemlab_module,'LinearElastostatics')
    if strcmp(config.vemlab_method,'VEM2D')
      [stresses,strains]=...
        vem_stresses_and_strains_linelast2d(u_nodal_sol,mesh,props,config);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      [stresses,strains]=...
        fem_stresses_and_strains_linelast2d(u_nodal_sol,mesh,props,config);     
    else
      throw_error('In compute_stresses_and_strains.m: vemlab_method\n');
    end
  else
    throw_error('In compute_stresses_and_strains.m: vemlab_module\n');
  end
end

