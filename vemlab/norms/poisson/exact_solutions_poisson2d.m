function [u_exact,flux_exact,grad_exact] = ...
              exact_solutions_poisson2d(exact_sol,eval_coords,matProps)

  % "eval_coords" can be  an integration point, a node, an array of integration 
  % points or an array of nodes 
  x=eval_coords(:,1);
  y=eval_coords(:,2);
  u_exact=exact_sol.u(x,y);  % [u1;u2;...;un]   
  grad_exact.x=exact_sol.dudx(x,y);
  grad_exact.y=exact_sol.dudy(x,y); 
  k=matProps.k;  % isotropic material  
  flux_exact.x=-k*exact_sol.dudx(x,y);
  flux_exact.y=-k*exact_sol.dudy(x,y);    
  
end





