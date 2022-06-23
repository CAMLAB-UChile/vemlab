function postprocess_numerical_solution_poisson2d(domainMesh,solution,matProps,config)
  % plot numerical solution
  [flux,grad]=...
      plot_numerical_solution_poisson2d(domainMesh,solution,matProps,config);                             
  % write numerical solutions to a text file
  if strcmp(config.write_solutions_to_text_file,'yes')
    write_solution_txt_poisson2d(domainMesh,solution,flux,grad,config);
  end
  % write numerical solutions to a CPP text file
  if strcmp(config.write_solutions_to_CPP_file,'yes')
    write_solution_CPP_poisson2d(domainMesh,solution,flux,grad,config);
  end  
  % write numerical solutions to a GiD file
  if strcmp(config.write_solutions_to_GiD_file,'yes')
    write_solution_GiD_poisson2d(domainMesh,solution,flux,grad,config);
  end
  % write numerical solutions to a VTK file
  if strcmp(config.write_solutions_to_VTK_file,'yes')
    write_solution_VTK_poisson2d(domainMesh,solution,config);
  end   
end

