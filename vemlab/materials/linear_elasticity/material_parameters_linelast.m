function matProps = material_parameters_linelast(Ey,nu,plane_state)
  %fprintf('Defining linear elastic material...\n'); 
  
  %
  % NOTE: THIS FUNCTION DEFINES D THAT IS CONSISTENT WITH A STRAIN VECTOR
  % DEFINED AS [e11; e22; e12].... NOTE THAT THE LAST TERM IS "e12"
  % WITHOUT 2 IN FRONT OF IT
  %
  
  matProps.Ey=Ey;  
  matProps.nu=nu;    
  matProps.plane_state=plane_state;
  if strcmp(plane_state,'plane_strain')
    lam=Ey*nu/((1+nu)*(1-2*nu));
    mu=Ey/(2*(1+nu));    
    matProps.D=(Ey/((1+nu)*(1-2*nu)))*[1-nu,nu,0; nu,1-nu,0; 0,0,2*(1-2*nu)];
    matProps.nu_bar=nu/(1-nu);
    matProps.Ey_bar=Ey/(1-nu*nu);   
    matProps.lam=lam;
    matProps.mu=mu;      
  elseif strcmp(plane_state,'plane_stress')
    lam=Ey*nu/((1+nu)*(1-2*nu));
    mu=Ey/(2*(1+nu));    
    matProps.D=(Ey/(1-nu*nu))*[1,nu,0; nu,1,0; 0,0,2*(1-nu)];
    matProps.nu_bar=nu;
    matProps.Ey_bar=Ey;    
    matProps.lam=lam;
    matProps.mu=mu;    
  else
    throw_error('Error in material_parameters.m: plane_state\n');    
  end    
end

