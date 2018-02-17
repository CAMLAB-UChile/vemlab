function [u_exact,dudx_exact,dudy_exact] = ...
              exact_solutions_linelast2d(exact_sol,eval_coords)

  % "eval_coords" can be  an integration point, a node, an array of integration 
  % points or an array of nodes
  num_nodes=length(eval_coords(:,1));
  u_exact=zeros(2*num_nodes,1);     % [u1x;u1y;u2x;u2y;...;unx;uny]
  dudx_exact=zeros(2*num_nodes,1);  % [du1x/dx;du1y/dx;du2x/dx;du2y/dx;...;dunx/dx;duny/dx]
  dudy_exact=zeros(2*num_nodes,1);  % [du1x/dy;du1y/dy;du2x/dy;du2y/dy;...;dunx/dy;duny/dy]  
  x=eval_coords(:,1);
  y=eval_coords(:,2);
  range=1:num_nodes;
  u_exact(2*range-1)=exact_sol.ux(x,y);    
  u_exact(2*range)=exact_sol.uy(x,y); 
  dudx_exact(2*range-1)=exact_sol.duxdx(x,y); 
  dudx_exact(2*range)=exact_sol.duydx(x,y);
  dudy_exact(2*range-1)=exact_sol.duxdy(x,y); 
  dudy_exact(2*range)=exact_sol.duydy(x,y);      

end





