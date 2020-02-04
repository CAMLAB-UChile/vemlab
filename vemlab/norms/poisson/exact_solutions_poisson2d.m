%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      VEMLab 
%-----------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Jan. 31, 2020: add a check on matProps to figure out if comes on an
%                element-by-element fashion or not.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u_exact,flux_exact,grad_exact] = ...
              exact_solutions_poisson2d(exact_sol,eval_coords,matProps)
            
  % figure out whether material's data is given in an element-by-element fashion or not
  size_k=length(matProps.k);
  if size_k > 1
    throw_warning('In exact_solutions_poisson2d.m : elements with possibly different conductivities.... assuming every element has the conductivity of the first one');  
    k=matProps.k{1,1}; % conductivity is particular for the current element  % isotropic material 
  else
    k=matProps.k;  % isotropic material
  end
  if length(k) > 1
    throw_error('In exact_solutions_poisson2d.m : multiple conductivities... not implemented for this condition');
  end  
  % "eval_coords" can be  an integration point, a node, an array of integration 
  % points or an array of nodes 
  x=eval_coords(:,1);
  y=eval_coords(:,2);
  u_exact=exact_sol.u(x,y);  % [u1;u2;...;un]   
  grad_exact.x=exact_sol.dudx(x,y);
  grad_exact.y=exact_sol.dudy(x,y);  
  flux_exact.x=-k*exact_sol.dudx(x,y);
  flux_exact.y=-k*exact_sol.dudy(x,y);    
  
end





