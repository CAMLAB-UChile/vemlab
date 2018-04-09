function write_solution_GiD_linelast2d(mesh,displacements,stresses,strains,config)

  if strcmp(config.vemlab_method,'VEM2D')||strcmp(config.vemlab_method,'FEM2DT3')
    write_solution_GiD_VEM2D_FEM2DT3_linelast2d(mesh,displacements,stresses,...
                                                strains,config);
  elseif strcmp(config.vemlab_method,'FEM2DQ4')
    write_solution_GiD_FEM2DQ4_linelast2d(mesh,displacements,stresses,...
                                          strains,config);    
  else
    throw_error('Error in write_solution_txt_linelast2d.m: vemlab_method')
  end

end

