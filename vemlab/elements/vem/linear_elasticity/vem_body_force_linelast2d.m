function fb=vem_body_force_linelast2d(verts,body_force_fun_values)
  % area of the element  
%   area_components=verts(:,1).*verts([2:end,1],2)-verts([2:end,1],1).*verts(:,2);
%   area=0.5*abs(sum(area_components));
  area=polyarea(verts(:,1),verts(:,2));  
  % number of vertices
  n=length(verts);
  % vector of nodal body forces
  bf=body_force_function_linelast2d(verts,body_force_fun_values);
  % mean body forces: this should be replaced by b_j = 1/|E| * int_E b_j dx
  mean_bfx=(1/n)*sum(bf(1:2:(2*n-1),1));
  mean_bfy=(1/n)*sum(bf(2:2:(2*n),1));
  mean_bf=[mean_bfx; mean_bfy];
  % VEM body force vector
  N_bar_i=[1/n,0;0,1/n];
  N_bar=repmat(N_bar_i,n,1);
  fb=area*N_bar*mean_bf;
  % this also works: here bf is a matrix containing the nodal values of bf
  % N_bar(1,2*n)=1/n;
  % fb=area*N_bar'.*bf;  
end