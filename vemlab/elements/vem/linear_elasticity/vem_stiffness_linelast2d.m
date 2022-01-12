function [Kvem]=vem_stiffness_linelast2d(verts,matProps,config)
  % area of the element
%   area_components=verts(:,1).*verts([2:end,1],2)-verts([2:end,1],1).*verts(:,2);
%   area=0.5*abs(sum(area_components));
  area=polyarea(verts(:,1),verts(:,2));
  % VEM element matrices
  Wc=vem_Wc_linelast2d(verts);
  Wr=vem_Wr_linelast2d(verts);
  Hc=vem_Hc_linelast2d(verts);
  Hr=vem_Hr_linelast2d(verts);
  Pp=Hr*Wr'+Hc*Wc';
  I_2N=eye(2*length(verts));
  
  nvertices=length(verts);
  if config.stability_type==1 % Gain et al.
    gamma=1.0;
    alpha=gamma*area*trace(matProps.D)/trace(Hc'*Hc);
    Se=alpha*I_2N;
  elseif config.stability_type==2  % D-recipe
    sizeKc=2*nvertices;    
    Se = max(eye(sizeKc),diag(diag(area*Wc*(matProps.D)*Wc')));   
  elseif config.stability_type==3  % modified D-recipe
    sizeKc=2*nvertices;
    Se=max((trace(Wc*(matProps.D)*(Wc')*area)/sizeKc)*eye(sizeKc),diag(diag(area*Wc*(matProps.D)*Wc')));   
  else
    throw_error('In vem_stiffness_linelast2d.m: stability_type');
  end

  if config.stability_type==0 % no stability
    Kvem=Wc*(matProps.D)*(Wc')*area; 
  else
    Kvem=Wc*(matProps.D)*(Wc')*area + (I_2N-Pp)'*Se*(I_2N-Pp); 
  end
  
end