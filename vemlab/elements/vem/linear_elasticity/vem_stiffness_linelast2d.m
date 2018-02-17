function [Kvem]=vem_stiffness_linelast2d(verts,matProps)
  % area of the element
  area_components=verts(:,1).*verts([2:end,1],2)-verts([2:end,1],1).*verts(:,2);
  area=0.5*abs(sum(area_components));
  % VEM element matrices
  Wc=vem_Wc_linelast2d(verts);
  Wr=vem_Wr_linelast2d(verts);
  Hc=vem_Hc_linelast2d(verts);
  Hr=vem_Hr_linelast2d(verts);
  Pp=Hr*Wr'+Hc*Wc';
  I_2N=eye(2*length(verts));
  % scaling parameter
  gamma=1.0;
  alpha=gamma*area*trace(matProps.D)/trace(Hc'*Hc);
  % VEM element stiffness
  Se=alpha*I_2N;
  Kvem=area*Wc*(matProps.D)*Wc'+(I_2N-Pp)'*Se*(I_2N-Pp);
end