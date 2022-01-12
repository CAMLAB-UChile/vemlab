function [u_exact,p_exact,strainvec_exact] = ...
              exact_solutions_linelast2d(exact_sol,eval_coords,matProps)

  % "eval_coords" can be  an integration point, a node, an array of integration 
  % points or an array of nodes
  num_nodes=length(eval_coords(:,1));
  u_exact=zeros(2*num_nodes,1);     % [u1x;u1y;u2x;u2y;...;unx;uny] 
  p_exact=zeros(num_nodes,1);       % [p1;p2;...;pn];  
  x=eval_coords(:,1);
  y=eval_coords(:,2);
  range=1:num_nodes;
  u_exact(2*range-1)=exact_sol.ux(x,y,matProps);    
  u_exact(2*range)=exact_sol.uy(x,y,matProps); 
  p_exact(range)=exact_sol.p(x,y,matProps);
  strainvec_exact=exact_sol.strainvec(x,y,matProps); % [e11_1 e22_1 e12_1;...;e11_n e22_n e12_n]   

end





