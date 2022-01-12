function bf = body_force_function_linelast2d(eval_coords,body_force_fun_values,matProps)
  % "eval_coords" can be  an integration point, a node, an array of integration 
  % points or an array of nodes
  num_nodes=size(eval_coords,1);
  x=eval_coords(:,1); 
  y=eval_coords(:,2);  
  range=1:num_nodes; 
  bf(2*range-1,1)=body_force_fun_values.bx(x,y,matProps); % bx
  bf(2*range,1)=body_force_fun_values.by(x,y,matProps); % by 
end

