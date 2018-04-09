function postprocess_numerical_solution_linelast2d(domainMesh,solution,matProps,config)
  % plot numerical solution
  [stress,strain]=plot_numerical_solution_linelast2d(domainMesh,solution,...
                                                     matProps,config);                             
  % write numerical solutions to a text file
  if strcmp(config.write_solutions_to_text_file,'yes')
    write_solution_txt_linelast2d(domainMesh,solution,stress,strain,config);
  end
  % write numerical solutions to a GiD file
  if strcmp(config.write_solutions_to_GiD_file,'yes')
    write_solution_GiD_linelast2d(domainMesh,solution,stress,strain,config);
  end
  % write numerical solutions to a VTK file
  if strcmp(config.write_solutions_to_VTK_file,'yes')
    write_solution_VTK_linelast2d(domainMesh,solution,config);
  end  
end

