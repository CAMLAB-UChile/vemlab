function [K_global,f_global] = vem_sparse_assembly_linelast2d(domainMesh,matProps,body_force_fun_values)
  num_nodes=length(domainMesh.coords(:,1)); 
  numel=length(domainMesh.connect);  
  
%   K_global=sparse(2*num_nodes,2*num_nodes);
%   f_global=sparse(2*num_nodes,1);
  
  % Lengths of buffer vectors for sparse assembly
  len_vector_K = 0; 
  len_vector_f = 0;
  % Determine exact lengths of buffer vectors
  % (Elements can have arbitrary numbers of nodes)
  for e=1:numel
      num_eldofs = 2*length(domainMesh.connect{e,1});
      len_vector_K = len_vector_K + (num_eldofs)^2; % unraveled size of K_local
      len_vector_f = len_vector_f + num_eldofs;
  end
  % Vectors for sparse assmbly of K_global
  vector_K = zeros(len_vector_K,1);
  indrow_K = zeros(len_vector_K,1);
  indcol_K = zeros(len_vector_K,1);
  % Vectors for sparse assmbly of f_global
  vector_f = zeros(len_vector_f,1);
  indrow_f = zeros(len_vector_f,1);
  pos_vector_K = 0;
  pos_vector_f = 0;
  for e=1:numel 
    nodes=domainMesh.connect{e,1};
    verts=domainMesh.coords(nodes,:);
    ind=zeros(2*length(nodes),1);      
    range=1:length(nodes);
    ind(2*range-1)=2*nodes-1;
    ind(2*range)=2*nodes;
    
%     K_global(ind,ind)=K_global(ind,ind)+vem_stiffness_linelast2d(verts,matProps);
%     f_global(ind)=f_global(ind)+vem_body_force_linelast2d(verts,body_force_fun_values);
         
    num_eldofs = 2*length(nodes);
    K_local = vem_stiffness_linelast2d(verts,matProps);
    f_local = vem_body_force_linelast2d(verts,body_force_fun_values);
    ind_vector_K = (1:num_eldofs^2) + pos_vector_K;
    ind_vector_f = (1:num_eldofs) + pos_vector_f;
    vector_K(ind_vector_K) = K_local(:);
    vector_f(ind_vector_f) = f_local;
    indrow_K_ = ind(:,ones(length(ind),1));
    indcol_K_ = indrow_K_';
    indrow_K(ind_vector_K) = indrow_K_(:);
    indcol_K(ind_vector_K) = indcol_K_(:);
    indrow_f(ind_vector_f) = ind;
    pos_vector_K = pos_vector_K + num_eldofs^2;
    pos_vector_f = pos_vector_f + num_eldofs;
  end
  K_global = sparse(indrow_K,indcol_K,vector_K,num_nodes*2,num_nodes*2);
  f_global = sparse(indrow_f,ones(len_vector_f,1),vector_f,num_nodes*2,1);
end

