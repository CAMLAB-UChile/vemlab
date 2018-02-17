function matProps = material_parameters_linelast(Ey,nu,plane_state)
  fprintf('Defining linear elastic material...\n'); 
  matProps.Ey=Ey;  
  matProps.nu=nu;    
  matProps.plane_state=plane_state;
  if strcmp(plane_state,'plane_strain')
    matProps.D=(Ey/((1+nu)*(1-2*nu)))*[1-nu,nu,0; nu,1-nu,0; 0,0,2*(1-2*nu)];
    matProps.nu_bar=nu/(1-nu);
    matProps.Ey_bar=Ey/(1-nu*nu);    
  elseif strcmp(plane_state,'plane_stress')
    matProps.D=(Ey/(1-nu*nu))*[1,nu,0; nu,1,0; 0,0,2*(1-nu)];
    matProps.nu_bar=nu;
    matProps.Ey_bar=Ey;    
  else
    throw_error('Error in material_parameters.m: plane_state\n');    
  end    
end

