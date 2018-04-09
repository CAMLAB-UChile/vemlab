function bf = source_term_function_poisson2d(eval_coords,source_term_fun_values)
  % "eval_coords" can be  an integration point, a node, an array of integration 
  % points or an array of nodes
  num_nodes=size(eval_coords,1);
  x=eval_coords(:,1); 
  y=eval_coords(:,2);  
  range=1:num_nodes; 
  bf(range,1)=source_term_fun_values(x,y); % bx
end

