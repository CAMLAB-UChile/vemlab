function [Kvem]=vem_stiffness_poisson2d(verts,matProps)
  % area of the element
  area_components=verts(:,1).*verts([2:end,1],2)-verts([2:end,1],1).*verts(:,2);
  area=0.5*abs(sum(area_components));
  % VEM element matrices
  Wc=vem_Wc_poisson2d(verts);
  Wr=vem_Wr_poisson2d(verts);
  Hc=vem_Hc_poisson2d(verts);
  Hr=vem_Hr_poisson2d(verts);
  Pp=Hr*Wr'+Hc*Wc';
  I_N=eye(length(verts));
  % VEM element stiffness
  k=matProps.k;    % isotropic material   
  Kvem=k*area*(Wc*Wc') + k*((I_N-Pp)')*(I_N-Pp); % not sure if k should multiply the stability part!
end