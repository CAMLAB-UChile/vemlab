function [h_max,L2rel,H1rel,L2prel]=compute_norms_of_the_error(exact_sol,uh_global,domainMesh,config,props)
  if strcmp(config.vemlab_module,'LinearElastostatics')
    if strcmp(config.vemlab_method,'VEM2D')
      [h_max,L2rel,H1rel,L2prel]=...
        vem_norms_linelast2d(exact_sol,uh_global,domainMesh,config,props);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      [h_max,L2rel,H1rel,L2prel]=...
        fem_norms_linelast2d(exact_sol,uh_global,domainMesh,config,props);      
    else
      throw_error('In compute_norms_of_the_error.m: vemlab_method\n');
    end
  elseif strcmp(config.vemlab_module,'Poisson')
    if strcmp(config.vemlab_method,'VEM2D')
      [h_max,L2rel,H1rel,L2prel]=vem_norms_poisson2d(exact_sol,uh_global,domainMesh,config,props);
    elseif strcmp(config.vemlab_method,'FEM2DT3')||strcmp(config.vemlab_method,'FEM2DQ4')
      [h_max,L2rel,H1rel,L2prel]=fem_norms_poisson2d(exact_sol,uh_global,domainMesh,config,props);      
    else
      throw_error('In compute_norms_of_the_error.m: vemlab_method\n');
    end    
  else
    throw_error('In compute_norms_of_the_error.m: vemlab_module\n');
  end
end

