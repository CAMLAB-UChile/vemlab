function fb=vem_source_vector_poisson2d(verts,source_term_fun_values)
  % area of the element  
  area_components=verts(:,1).*verts([2:end,1],2)-verts([2:end,1],1).*verts(:,2);
  area=0.5*abs(sum(area_components));
  % number of vertices
  N=length(verts);
  % vector of nodal source terms
  b=source_term_function_poisson2d(verts,source_term_fun_values);
  % mean body forces: this should be replaced by b_j = 1/|E| * int_E b_j dx
  mean_b=(1/N)*sum(b);
  % VEM source vector
  N_bar=(1/N)*ones(N,1);
  fb=area*N_bar*mean_b;
end